#include "FE_sweeper.hpp"

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
            
        std::shared_ptr<VectorType> w; 

            
        public:
            explicit fischer_sweeper<SweeperTrait, BaseFunction, Enabled>(std::shared_ptr<BaseFunction> basis, size_t nlevel, std::shared_ptr<GridType> grid)
                                    : Heat_FE<SweeperTrait, BaseFunction, Enabled>(basis, nlevel, grid){
        
                this->nlevel = nlevel;    
	    
                this->basis = basis;

                this->grid = grid;
	
                assembleProblem(basis, this->A_dune, this->M_dune);

                this->stiffnessMatrix = this->A_dune;
                this->stiffnessMatrix *= -1;
        
                w = std::make_shared<VectorType>(this->M_dune.M());
        
                for(int j=0; j<this->M_dune.M(); j++){(*w)[j]=0;}

                for(int i=0; i<this->M_dune.M(); i++){
                    for(int j=0; j<this->M_dune.M(); j++){
                        if(this->M_dune.exists(i,j))
                            (*w)[i][0]= ((double) (*w)[i][0]) + ((double) this->M_dune[i][j][0][0]);
                    }
                }

                const auto bs = basis->size();
                std::cout << "Finite Element basis of level " << nlevel << " consists of " <<  basis->size() << " elements " << std::endl;

                this->encap_factory()->set_size(bs);

          }
            
            
            
        
          fischer_sweeper(const fischer_sweeper<SweeperTrait, BaseFunction, Enabled>& other) = default;
          fischer_sweeper(fischer_sweeper<SweeperTrait, BaseFunction, Enabled>&& other) = default;
          virtual ~fischer_sweeper() = default;
            
          shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t)      {
            auto result = this->get_encap_factory().create();
            const auto dim = 1; //SweeperTrait::DIM;
            double n  = this-> _n;
            double l0 = this-> _nu;
            double l1 = l0/2. *(pow((1+n/2.), 1/2.) + pow((1+ n/2.), -1/2.) );
            double d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
            auto exact_solution = [l0, l1, n, d, t](const Dune::FieldVector<double,dim>&x){ 
                return pow((1 + (pow(2, n/2.)-1 )* exp(-(n/2.)*d*(x+2*l1*t)) ), -2./n);
            };  
            auto N_x = [t](const Dune::FieldVector<double,dim>&x){
                return x;
            };
            Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
            interpolate(*this->basis, x_node, N_x);
            interpolate(*this->basis, result->data(), exact_solution);
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
