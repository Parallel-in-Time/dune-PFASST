#include "FE_sweeper.hpp"
#include "functions.hh"


#include <dune/solvers/iterationsteps/blockgssteps.hh>
#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/common/defaultbitvector.hh>
#include <dune/solvers/common/resize.hh>
#include <dune/solvers/transferoperators/compressedmultigridtransfer.hh>
#include <dune/solvers/transferoperators/densemultigridtransfer.hh>
#include <dune/solvers/iterationsteps/multigridstep.hh>
#include <dune/solvers/solvers/umfpacksolver.hh>
#include <dune/solvers/norms/energynorm.hh>
#include <dune/solvers/norms/twonorm.hh>


#include <dune/solvers/iterationsteps/blockgssteps.hh>
#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/common/defaultbitvector.hh>
#include <dune/solvers/common/resize.hh>
#include <dune/solvers/transferoperators/compressedmultigridtransfer.hh>
#include <dune/solvers/iterationsteps/multigridstep.hh>
#include <dune/solvers/solvers/umfpacksolver.hh>
#include <dune/solvers/norms/energynorm.hh>
#include <dune/solvers/norms/twonorm.hh>


#include <dune/tnnmg/iterationsteps/tnnmgstep.hh>
#include <dune/tnnmg/iterationsteps/nonlineargsstep.hh>
#include <dune/tnnmg/functionals/boxconstrainedquadraticfunctional.hh>
#include <dune/tnnmg/functionals/bcqfconstrainedlinearization.hh>
#include <dune/tnnmg/projections/obstacledefectprojection.hh>
#include <dune/tnnmg/localsolvers/scalarobstaclesolver.hh>

//#include <dune/tnnmg/localsolvers.hh>

// eigene TNNMG Funktionale
#include "tnnmgfunctional.hh"
#include "tnnmgfunctionallinearization.hh"
#include "scalarbisectionsolver.hh"


#include<iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/schwarz.hh>
//#include <dune/istl/matrixmarket.hh>
#include<dune/istl/paamg/pinfo.hh>
//#include<dune/istl/matrixredistribute.hh>
#include<dune/istl/paamg/graph.hh>

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
        double                                     	_nu{0.0};
        double                                     	_n{1.0};
        double                                      	_delta{1.0};
        double                                          _abs_newton_tol=1e-10; 
            
        public:
            explicit fischer_sweeper<SweeperTrait, BaseFunction, Enabled>(std::shared_ptr<BaseFunction> basis, size_t nlevel, std::shared_ptr<GridType> grid)
                                    : Heat_FE<SweeperTrait, BaseFunction, Enabled>(basis, nlevel, grid){
        
                //using MatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double, blockSize, blockSize> >;
                //MatrixType mass;
                //MatrixType stiffness;
                std::cout << "***************************************************** Konstruktor anfang *********************************************************************************" << std::endl;

                    MatrixType ipdg;
                    auto iBackend = Dune::Fufem::istlMatrixBackend(ipdg);
                    auto sBackend = Dune::Fufem::istlMatrixBackend(this->A_dune);
                    auto mBackend = Dune::Fufem::istlMatrixBackend(this->M_dune);

                    using Assembler = Dune::Fufem::DuneFunctionsOperatorAssembler<BaseFunction, BaseFunction>;
                    auto assembler = Assembler{*(this->basis), *(this->basis)};

                    using FiniteElement = std::decay_t<decltype(this->basis->localView().tree().finiteElement())>;
 
                    auto vintageLaplace = LaplaceAssembler<GridType,FiniteElement, FiniteElement>();
                    auto vintageMass = MassAssembler<GridType,FiniteElement, FiniteElement>();
 
                    auto localAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                     vintageLaplace.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
                    };
 
                    // IPDG terms ...
                    auto vintageIPDGAssembler = InteriorPenaltyDGAssembler<GridType, FiniteElement, FiniteElement>();
                    vintageIPDGAssembler.sigma0=penalty;
                    vintageIPDGAssembler.dirichlet = false;
                    auto localIPDGAssembler = [&](const auto& edge, auto& matrixContainer,
                         auto&& insideTrialLocalView, auto&& insideAnsatzLocalView, auto&& outsideTrialLocalView, auto&& outsideAnsatzLocalView)
                     {
                         vintageIPDGAssembler.assembleBlockwise(edge, matrixContainer, insideTrialLocalView.tree().finiteElement(),
                                                             insideAnsatzLocalView.tree().finiteElement(),
                                                             outsideTrialLocalView.tree().finiteElement(),
                                                             outsideAnsatzLocalView.tree().finiteElement());
                     };
                     auto localBoundaryAssembler = [&](const auto& edge, auto& localMatrix, auto&& insideTrialLocalView, auto&& insideAnsatzLocalView)
                     {
                     vintageIPDGAssembler.assemble(edge, localMatrix, insideTrialLocalView.tree().finiteElement(), insideAnsatzLocalView.tree().finiteElement());
                     };
 
 
                     auto patternBuilder = iBackend.patternBuilder();
                     assembler.assembleSkeletonPattern(patternBuilder);
                     patternBuilder.setupMatrix();
 
                     auto localMassAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                     vintageMass.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
                     };
 
                     assembler.assembleBulk(sBackend, localAssembler);
                     assembler.assembleBulk(mBackend, localMassAssembler);
                     assembler.assembleSkeletonEntries(iBackend, localIPDGAssembler, localBoundaryAssembler); // assemble IPDG terms into ipdg
                     ipdg+=this->A_dune;
                     this->A_dune = ipdg; // fix this TODO
                
                                        
                                        
                w = std::make_shared<VectorType>(this->basis->size());//this->M_dune.M());
        
                VectorType tmp(this->basis->size()); 
                tmp=1;
                this->M_dune.mv(tmp, *w); // w is just the row sum of mass
                std::cout << "***************************************************** Konstruktor ende *********************************************************************************" << std::endl;

          }
            
            
            
        
          fischer_sweeper(const fischer_sweeper<SweeperTrait, BaseFunction, Enabled>& other) = default;
          fischer_sweeper(fischer_sweeper<SweeperTrait, BaseFunction, Enabled>&& other) = default;
          virtual ~fischer_sweeper() = default;
            
          shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t)      {
          auto result = this->get_encap_factory().create();
            
          const int dim=2;
          ScalarRandomDiscFunction<dim> initialPhaseField(29007122, 50, 0.1, 0.5, -1.0, 1.0);
          auto phaseField = [&](auto x) {Dune::FieldVector<double, 1> y;initialPhaseField.evaluate(x,y); return y;};
          Dune::Functions::interpolate(*(this->basis), Dune::Fufem::istlVectorBackend(result->data()), phaseField);
             
//             auto result = this->get_encap_factory().create();
//             const auto dim = 1; //SweeperTrait::DIM;
//             double n  = this-> _n;
//             double l0 = this-> _nu;
//             double l1 = l0/2. *(pow((1+n/2.), 1/2.) + pow((1+ n/2.), -1/2.) );
//             double d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
//             auto exact_solution = [l0, l1, n, d, t](const Dune::FieldVector<double,dim>&x){ 
//                 return pow((1 + (pow(2, n/2.)-1 )* exp(-(n/2.)*d*(x+2*l1*t)) ), -2./n);
//             };  
//             interpolate(*this->basis, result->data(), exact_solution);
            
            return result;
          }
          
          shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,const shared_ptr<typename SweeperTrait::encap_t> u) {
	
            ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);
            auto result = this->get_encap_factory().create();
            double nu =this->_nu;
            for (int i=0; i<u->get_data().size(); ++i)
                for(int j=0; j<blockSize; j++){
                result->data()[i][j]= -pow(u->get_data()[i][j], this->_n+1) * (*w)[i][j];
                //std::cout << result->data()[i][j] << std::endl;
                }

            this->M_dune.umv(u->get_data(), result->data());
            result->data()*=nu*nu;
            this->A_dune.mmv(u->get_data(), result->data()); //umv

          	  std::cout << "evaluate "  << std::endl;  
//             for (size_t i = 0; i < u->get_data().size(); i++) {
//                 for(int j=0; j< blockSize; j++)
//                 std::cout << "f evaluate" << result->data()[i][j] << std::endl;
//             }
            return result;
        }
        
        
        
        void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs){
                                                        

            std::cout << "impl solve" << std::endl;  
            auto rhs_copy = this->get_encap_factory().create();
            rhs_copy = rhs;
            std::cout << "impl solve rhs " << rhs->data()[0][0] << std::endl;  
            
            ML_CVLOG(4, this->get_logger_id(),"IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);

            auto phiprime =
             [&](auto&& uu) {
             return dt*_nu*_nu*pow(uu, _n+1);
            }; 
            auto phi2prime = 
             [&] (auto&& uu) {
             return dt*(_nu*_nu)*(_n+1) * pow(uu, _n);
            };
 	         
            auto matrix_ = this->A_dune; // avoid too small matrix pattern
            matrix_=0;
            matrix_+=this->M_dune;
            matrix_ *= (1-dt*_nu*_nu);
            matrix_.axpy(dt, this->A_dune); 
   
            //std::cout << "Anfang TNNMG" << std::endl;
            using BitVector = Dune::Solvers::DefaultBitVector_t<VectorType>;
            BitVector ignore(u->data().size());
 
            std::cout << "vor tranfer" << std::endl;  
            using TransferOperator = DenseMultigridTransfer<VectorType>;
            using TransferOperators = std::vector<std::shared_ptr<TransferOperator>>;
            using TransferMatrices = std::vector<std::shared_ptr<MatrixType>>;

            // Set up dg transfer
            TransferMatrices transferMatrices;
            auto& gridptr = this->grid;         
            transferMatrices.resize(gridptr->maxLevel());
            for (auto& tt: transferMatrices)
                tt = std::make_shared<MatrixType>();
            Dune::HPDG::assembleDGGridTransferHierarchy(transferMatrices, *(this->grid));
            TransferOperators transfer(this->grid->maxLevel());
            for (size_t i = 0; i < transfer.size(); ++i)
            {
                // create transfer operator from level i to i+1
                transfer[i] = std::make_shared<TransferOperator>();
                std::cout << "Transfer size for level " << i <<": " << transferMatrices[i]->N() << "x" << transferMatrices[i]->M() << std::endl;
                transfer[i]->setMatrix(*transferMatrices[i]);
            }
            //// TNNMG without actual constraints
            using Functional = Dune::TNNMG::EnergyFunctional<MatrixType, VectorType, decltype(phiprime), decltype(phiprime), decltype(phi2prime), double>;
            auto J = Functional(matrix_, rhs->data(), *(this->w), phiprime, phiprime, phi2prime);
            auto localSolver = gaussSeidelLocalSolver(Dune::TNNMG::ScalarBisectionSolver());
            using NonlinearSmoother = Dune::TNNMG::NonlinearGSStep<Functional, decltype(localSolver), BitVector>;
            auto nonlinearSmoother = std::make_shared<NonlinearSmoother>(J, u->data(), localSolver);
            using Linearization = Dune::TNNMG::EnergyFunctionalConstrainedLinearization<Functional, BitVector>; // incorporates what used to be evaluate_f, evaluate_df
            using LineSearchSolver = Dune::TNNMG::ScalarBisectionSolver;
            /* Setup linear multigrid */
            using MultiGrid =Dune::Solvers::MultigridStep<MatrixType, VectorType, BitVector>;
            auto mgStep = std::make_shared<MultiGrid>();
            auto gssmoother = Dune::Solvers::BlockGSStepFactory<MatrixType, VectorType, BitVector>::create(Dune::Solvers::BlockGS::LocalSolvers::gs());
            mgStep->setSmoother(&gssmoother);
            mgStep->setTransferOperators(transfer);
            mgStep->setMGType(1,3,3);

            // base solver for multigrid
            auto umfpack = Dune::Solvers::UMFPackSolver<MatrixType, VectorType>{}; // direct solver
            auto trivialProjection = [](auto&& f, auto& x, auto& c) {};
            using Step = Dune::TNNMG::TNNMGStep<Functional, BitVector, Linearization, decltype(trivialProjection), LineSearchSolver>;
            int mu=1;
            // J Functional to minimize, mgStep: linear "solver" for linear correction, mu #multigrid steps per iteration, trivialProjection is identity
            auto step = Step(J, u->data(), nonlinearSmoother, mgStep, mu, trivialProjection, LineSearchSolver());
            step.setPreSmoothingSteps(0);
            step.setIgnore(ignore);

            //using Norm =  TwoNorm<VectorType>;
            using Norm =  EnergyNorm<MatrixType,VectorType>;
            auto norm = Norm(this->A_dune);
            // max. 20 iterations or two norm of correction less than 1e-10
            using Solver = LoopSolver<VectorType>;
            auto solver = Solver(&step, 20, 1e-10, &norm, Solver::FULL);

            solver.addCriterion(
                [&](){
                return Dune::formatString("   % 12.5e", step.lastDampingFactor());
                },
                "   damping     ");

            solver.addCriterion(
                [&](){
                return Dune::formatString("   % 12d", step.linearization().truncated().count());
                },
                "   truncated   ");

            solver.preprocess();
            solver.solve();

            Dune::BlockVector<Dune::FieldVector<double,blockSize> > M_u;
            M_u.resize(u->get_data().size());
            this->M_dune.mv(u->get_data(), M_u);
            //std::cout << "impl solve rhs_copy " << rhs_copy->data()[0][0] << std::endl;  

            //std::cout << "impl solve "  << std::endl;
            for (size_t i = 0; i < u->get_data().size(); i++) {
                for(int j=0; j<blockSize; j++){
                    f->data()[i][j] = (M_u[i][j] - rhs_copy->get_data()[i][j]) / (dt);
                    //std::cout << "f impl_solve" << f->data()[i][j] << std::endl;
                }
            }
            //evaluate_rhs_impl(0, u);
            //std::exit(0);
            this->_num_impl_solves++;
             
        }
        
        
        
        
        
        /*void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs){

         std::cout << "im implicit solve " <<  std::endl;            
   
         ML_CVLOG(4, this->get_logger_id(), "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);

         auto residuum = this->get_encap_factory().create();
	 Dune::BlockVector<Dune::FieldVector<double,blockSize> > newton_rhs, newton_rhs2 ;
         newton_rhs.resize(rhs->get_data().size());
         newton_rhs2.resize(rhs->get_data().size());
    
         u->zero();
	 for (int i=0; i< 200 ;i++){
            Dune::BCRSMatrix<Dune::FieldMatrix<double,blockSize,blockSize> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,blockSize,blockSize> >(this->A_dune); ///////M
            
            evaluate_f(f, u, dt, rhs);
            
            evaluate_df(df, u, dt);

            df.mv(u->data(), newton_rhs);
            newton_rhs -= f->data();
            newton_rhs2 = newton_rhs;
          
//             for (size_t i = 0; i < df.M(); i++) {
//             for (size_t j = 0; j < df.N(); j++)                 
//                 std::cout << "f " << df[i][j][0][0] << std::endl;
//             }
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
          for(int j=0; j< blockSize; j++){
          f->data()[i][j] = (M_u[i][j] - rhs->get_data()[i][j]) / (dt);
	  std::cout << "f implicit solve" << f->data()[i][j] << std::endl;}
        }
        evaluate_rhs_impl(0, u);
	std::exit(0);
        this->_num_impl_solves++;
      }*/
      
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
            for(int j=0; j< blockSize; j++){  
            f->data()[i][j]= pow(u->get_data()[i][j], _n+1) * (*w)[i][j];}	
          }
          this->M_dune.mmv(u->get_data(), f->data());

          f->data() *= (_nu*_nu);


          this->A_dune.mmv(u->get_data(),f->data());
          f->data() *= dt;
          this->M_dune.umv(u->get_data(),f->data());
          f->data() -=rhs->get_data();
      }
						
      void evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,blockSize,blockSize> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						){
            double _nu=this->_nu;
            double _n=this->_n;


            
            for (int i=0; i<u->get_data().size(); ++i)
            {
                for(int j=0; j< blockSize; j++)
                df[i][i][j][j]= (_nu*_nu)*(_n+1) * pow(u->get_data()[i][j], _n) * ((double) (*w)[i][j]);	
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
