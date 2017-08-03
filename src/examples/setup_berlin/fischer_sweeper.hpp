#include "FE_sweeper.hpp"



#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/istlbackend.hh>
#include <dune/fufem/formatstring.hh>




#include "functions.hh"

//this is a specified sweeper for the generalized Fisher equation
//it inherits from the implicit Finite Element sweeper

//you just need to add a few things

//the konstructer: mass- and stiffnessmatrix are already assembeld in the base class, but here we make a another discretisation vector w, because we use lumping for the nonlinearity

//exact: gives the exact solution of Fischers equation, we use it in the main to construct the initial value and to calculate the error after the simulation

//evaluate_rhs_impl: it gives back the right hand side of the discrtised ODE

//implicit_solve: it solves the resulting algebraic system which results from using sdc 

using namespace pfasst::examples::FE_sweeper;


namespace pfasst
{
  namespace examples
  {
    namespace fischer_example
    {
      template<
        class SweeperTrait,
        class BaseFunction,
        typename Enabled = void
      >
      class fischer_sweeper
        : public Heat_FE<SweeperTrait, BaseFunction, Enabled>{
            
        std::shared_ptr<VectorType>                     w; 
        double                                     	_nu{10};
        double                                     	_n{2.0};
        bool                                            _write = true;
        //double                                      	_delta{1.0};
        //double                                          _abs_newton_tol=1e-10; 
            
        public:
            explicit fischer_sweeper<SweeperTrait, BaseFunction, Enabled>(std::shared_ptr<BaseFunction> basis, size_t nlevel, std::shared_ptr<GridType> grid)
                                    : Heat_FE<SweeperTrait, BaseFunction, Enabled>(basis, nlevel, grid){
        
                auto sBackend = Dune::Fufem::istlMatrixBackend(this->A_dune);
                auto mBackend = Dune::Fufem::istlMatrixBackend(this->M_dune);

                using Assembler = Dune::Fufem::DuneFunctionsOperatorAssembler<BaseFunction, BaseFunction>;
                auto assembler = Assembler{*basis, *basis};

                using FiniteElement = std::decay_t<decltype(basis->localView().tree().finiteElement())>;

                auto vintageLaplace = LaplaceAssembler<GridType,FiniteElement, FiniteElement>();
                auto vintageMass = MassAssembler<GridType,FiniteElement, FiniteElement>();

                auto localAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageLaplace.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
                };
                auto localMassAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageMass.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
                };

                assembler.assembleBulk(sBackend, localAssembler);
                assembler.assembleBulk(mBackend, localMassAssembler);
                
                this->A_dune*=-1;
    
                w = std::make_shared<VectorType>(this->M_dune.M());
        
                for(int j=0; j<this->M_dune.M(); j++){(*w)[j]=0;}

                for(int i=0; i<this->M_dune.M(); i++){
                    for(int j=0; j<this->M_dune.M(); j++){
                        if(this->M_dune.exists(i,j))
                            (*w)[i][0]= ((double) (*w)[i][0]) + ((double) this->M_dune[i][j][0][0]);
                    }
                }


          }
            
            
            
        
          fischer_sweeper(const fischer_sweeper<SweeperTrait, BaseFunction, Enabled>& other) = default;
          fischer_sweeper(fischer_sweeper<SweeperTrait, BaseFunction, Enabled>&& other) = default;
          virtual ~fischer_sweeper() = default;
          
             //write_results
          void write_results(shared_ptr<typename SweeperTrait::encap_t> x, const std::string filename){
             //auto xBE = Dune::Fufem::istlVectorBackend(x->data());
             //auto xFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(this->basis, Dune::TypeTree::hybridTreePath(), xBE);
             //this->grid.leafGridView();
            // Dune::VTKWriter<typename GridType::LeafGridView> vtkWriter(this->grid.leafGridView());
// 
//             vtkWriter.addVertexData(xFunction, Dune::VTK::FieldInfo("x", Dune::VTK::FieldInfo::Type::scalar, 1));
//             vtkWriter.write(filename);
          }
          
          shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t)      {
            auto result = this->get_encap_factory().create();
            
             const int dim=2;
             ScalarRandomDiscFunction<dim> initialPhaseField(29007122, 50, 0.1, 0.5, -1.0, 1.0);
             auto phaseField = [&](auto x) {Dune::FieldVector<double, 1> y;initialPhaseField.evaluate(x,y); return y;};
             Dune::Functions::interpolate(*(this->basis), Dune::Fufem::istlVectorBackend(result->data()), phaseField);

//              const auto dim = 1; //SweeperTrait::DIM;
//              double n  = this-> _n;
//              double l0 = this-> _nu;
//              double l1 = l0/2. *(pow((1+n/2.), 1/2.) + pow((1+ n/2.), -1/2.) );
//              double d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
//              auto exact_solution = [l0, l1, n, d, t](const Dune::FieldVector<double,dim>&x){ 
//                  return pow((1 + (pow(2, n/2.)-1 )* exp(-(n/2.)*d*(x+2*l1*t)) ), -2./n);
//              };  
//              auto N_x = [t](const Dune::FieldVector<double,dim>&x){
//                  return x;
//              };
//              Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
//              interpolate(*this->basis, x_node, N_x);
//              interpolate(*this->basis, result->data(), exact_solution);
             return result;
           }
           
           shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,const shared_ptr<typename SweeperTrait::encap_t> u) {
 	
             ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);
             auto result = this->get_encap_factory().create();
             double nu =this->_nu;
             for (int i=0; i<u->get_data().size(); ++i)
                 result->data()[i]= -pow(u->get_data()[i], this->_n+1) * (*w)[i];	
             this->M_dune.umv(u->get_data(), result->data());
             result->data()*=nu*nu;
             this->A_dune.umv(u->get_data(), result->data());
        
            return result;
        }
        
        void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs){
	
         ML_CVLOG(4, this->get_logger_id(), "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);

         //for(int i=0; i< u->data().size(); i++) std::cout << u->data()[i] << std::endl;

         
         auto residuum = this->get_encap_factory().create();
	 Dune::BlockVector<Dune::FieldVector<double,1> > newton_rhs, newton_rhs2 ;
         newton_rhs.resize(rhs->get_data().size());
         newton_rhs2.resize(rhs->get_data().size());
    
         u->zero();
	 for (int i=0; i< 200 ;i++){
            Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune); ///////M
            evaluate_f(f, u, dt, rhs);
            evaluate_df(df, u, dt);
            df.mv(u->data(), newton_rhs);
            newton_rhs -= f->data();
            newton_rhs2 = newton_rhs;
          
            Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df);
	  
            Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df,1.0);
	  
            Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-16, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              1);    // verbosity of the solver
          
          
            Dune::InverseOperatorResult statistics ;
            cg.apply(u->data(), newton_rhs , statistics ); //rhs ist nicht constant!!!!!!!!!

            evaluate_f(f, u, dt, rhs);
          
            std::cout << i << " residuumsnorm von f(u) " << f->norm0() << std::endl;  
            if(f->norm0()<1e-10){   std::cout << "genauigkeit erreicht " << i << std::endl;      break;} //  std::exit(0); std::cout << "genauigkeit erreicht " << i << std::endl;
          
            df.mv(u->data(), residuum->data());
            residuum->data() -= newton_rhs2;
            std::cout << "residuums norm " << residuum->norm0() << std::endl;

	}

	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);

        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
	  //std::cout << "f " << f->data()[i] << std::endl;
        }
        //evaluate_rhs_impl(0, u);
	//std::exit(0);
        this->_num_impl_solves++;
      }
      
        private:
            
      void evaluate_f(shared_ptr<typename SweeperTrait::encap_t> f,
            const shared_ptr<typename SweeperTrait::encap_t> u,
            const typename SweeperTrait::time_t& dt,
            const shared_ptr<typename SweeperTrait::encap_t> rhs){

          double _nu=this->_nu;
          double _n=this->_n;

          f->zero();
          Dune::BlockVector<Dune::FieldVector<double,1> > fneu;
          fneu.resize(u->get_data().size());
          for (int i=0; i<u->get_data().size(); ++i)
          {
            f->data()[i]= pow(u->get_data()[i], _n+1) * (*w)[i];	
          }
          this->M_dune.mmv(u->get_data(), f->data());

          f->data() *= (_nu*_nu);


          this->A_dune.mmv(u->get_data(),f->data());
          f->data() *= dt;
          this->M_dune.umv(u->get_data(),f->data());
          f->data() -=rhs->get_data();
      }
						
      void evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						){
            double _nu=this->_nu;
            double _n=this->_n;
            for (int i=0; i<df.N(); ++i)
            {
                df[i][i]= (_nu*_nu)*(_n+1) * pow(u->get_data()[i], _n) * ((double) (*w)[i]);	
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
      } 
          
    };   
    }
  }
}
