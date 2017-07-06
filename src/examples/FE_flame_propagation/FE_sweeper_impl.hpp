#include "FE_sweeper.hpp"

#include "assemble.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <utility>
#include <vector>
using std::shared_ptr;
using std::vector;

//#include <leathers/push>
//#include <leathers/all>
//#include <boost/math/constants/constants.hpp>
//#include <leathers/pop>
//using boost::math::constants::pi;
//using boost::math::constants::two_pi;
//using boost::math::constants::pi_sqr;

double pi = 3.14159265359;
double two_pi = 2*pi;
double pi_sqr= pi*pi;


#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>

#include <iostream>

// dune solvers
#include <dune/solvers/iterationsteps/blockgssteps.hh>
#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/common/defaultbitvector.hh>
#include <dune/solvers/common/resize.hh>
#include <dune/solvers/transferoperators/compressedmultigridtransfer.hh>
#include <dune/solvers/iterationsteps/multigridstep.hh>
#include <dune/solvers/solvers/umfpacksolver.hh>

#include <dune/tnnmg/iterationsteps/tnnmgstep.hh>
#include <dune/tnnmg/iterationsteps/nonlineargsstep.hh>
#include <dune/tnnmg/functionals/boxconstrainedquadraticfunctional.hh>
#include <dune/tnnmg/functionals/bcqfconstrainedlinearization.hh>
#include <dune/tnnmg/projections/obstacledefectprojection.hh>
#include <dune/tnnmg/localsolvers/scalarobstaclesolver.hh>
//#include <dune/tnnmg/localsolvers.hh>
#include <dune/solvers/norms/energynorm.hh>


//#include <c++/4.8/memory>
struct TrivialSolver {
  template<class Vector, class Functional, class BitVector>
    constexpr void operator()(Vector& x, const Functional& f, const BitVector& ignore) const
    { x=1.0;}
};

namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::init_opts()
      {
        /*config::options::add_option<size_t>("Heat FE", "num_dofs",
                                            "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat FE", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat FE", "nu",
                                               "thermal diffusivity");*/
      }


        template<class SweeperTrait, typename Enabled>
      Heat_FE<SweeperTrait, Enabled>::Heat_FE(std::shared_ptr<fe_manager> FinEl, size_t nlevel)
        :   IMEX<SweeperTrait, Enabled>()

      {
      
	/*Dune::FieldVector<double,SweeperTrait::DIM> hR = {20};
	Dune::FieldVector<double,SweeperTrait::DIM> hL = {-20};
        array<int,SweeperTrait::DIM> n;

	std::fill(n.begin(), n.end(), nelements);

        this->grid = std::make_shared<GridType>(hL, hR, n);
        //grid.globalRefine(0);
	
        typedef GridType::LeafGridView GridView;
        GridType::LeafGridView gridView = grid->leafGridView();

        std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
        std::cout << "***** Ordnung der Basis " << SweeperTrait::BASE_ORDER << std::endl;
	
	
        this->basis = std::make_shared<BasisFunction>(gridView);

        std::cout << "***** Basis erstellt mit " <<  basis->size() << " Elementen " << std::endl;

        this->encap_factory()->set_size(basis->size());*/
	this->FinEl = FinEl;
	basis = FinEl->get_basis(nlevel);
	    
	assembleProblem(basis, this->A_dune, this->M_dune);

        const auto bs = basis->size();
        std::cout << "Finite Element basis of level " << nlevel << " consists of " <<  basis->size() << " elements " << std::endl;

        this->encap_factory()->set_size(bs);

      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::set_options()
      {


        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        int num_nodes = this->get_quadrature()->get_num_nodes();

        //assembleProblem(basis, A_dune, M_dune);


      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        const auto dim = 1; //SweeperTrait::DIM;
        spatial_t nu = this-> _nu; 
	spatial_t delta = this->_delta;
	//std::cout << "nu = " << this->_nu << std::endl;
	//std::cout << "delta = " << this->_delta << std::endl;
        auto exact_solution = [t, nu, dim, delta](const Dune::FieldVector<double,dim>&x){
          double c = 2./delta;  
	  return 0.5*(1.0-std::tanh((x[0] - c*t*nu)/(delta)));
        };


	

        auto N_x = [t](const Dune::FieldVector<double,dim>&x){
            return x;

        };

        Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
        interpolate(*basis, x_node, N_x);

        interpolate(*basis, result->data(), exact_solution);

	for (int i=0; i< result->data().size(); i++){
	//std::cout << "result = " << result->data()[i] << std::endl;
	}
        return result;
      }
      
      
      
      
      //______________________________
      
      
      
      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::source(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        const auto dim = 1;
        spatial_t nu = this-> _nu; // nu je valjda desni rub
	//std::cout << "nu = " << this->_nu << std::endl;
        auto exact_solution_source = [t, nu, dim](const Dune::FieldVector<double,dim>&x){
            double solution=1.0;
            //for(int i=0; i<SweeperTrait::DIM; i++){solution *=x[i];}    //
            
	    std::cout << "PI = " << PI << std::endl;
	    std::exit(0);
	    return 2*t + PI*PI*(std::sin(PI*x[0])*x[0] + std::sin(PI*x[1])*x[1])  - 2*PI*(std::cos(PI*x[0]) + std::cos(PI*x[1])) + pow(std::sin(PI*x[0])*x[0] + std::sin(PI*x[1])*x[1] + t*t, 2) ;
        };


	
	

        auto N_x = [t](const Dune::FieldVector<double,dim>&x){
            return x;

        };

        Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
        interpolate(*basis, x_node, N_x);

        interpolate(*basis, result->data(), exact_solution_source);


        return result;
      }
      
      
      
      
      
      
      
      
      
      //_______________________________
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, Enabled>::converged(pre_check);

        if (!pre_check) {
          assert(this->get_status() != nullptr);
          const typename traits::time_t t = this->get_status()->get_time();
          const typename traits::time_t dt = this->get_status()->get_dt();

          assert(this->get_quadrature() != nullptr);
          auto nodes = this->get_quadrature()->get_nodes();
          const auto num_nodes = this->get_quadrature()->get_num_nodes();
          nodes.insert(nodes.begin(), typename traits::time_t(0.0));

          ML_CVLOG(1, this->get_logger_id(),
                   "Observables after "
                   << ((this->get_status()->get_iteration() == 0)
                          ? std::string("prediction")
                          : std::string("iteration ") + std::to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(),
                     "  t["<<m<<"]=" <<  (t + dt * nodes[m])
                     << "      |abs residual| = " <<  this->_abs_res_norms[m]
                     << "      |rel residual| = " <<  this->_rel_res_norms[m]
//                      << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
//                      << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" <<  (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " <<  this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " <<  this->_rel_res_norms[num_nodes]
//                   << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
//                   << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes])
                 );
        }
        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      Heat_FE<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }

      typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; //ruth_dim
     
     
      
      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<GridType>
      Heat_FE<SweeperTrait, Enabled>::get_grid() const
      {
        return grid;
      }*/


      template<class SweeperTrait, typename Enabled> //Fehler der aktuellen Loesung an jedem Quadraturpunkt
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
      {
        ML_CVLOG(4, this->get_logger_id(), "computing error");

        assert(this->get_status() != nullptr);
        const typename traits::time_t dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), typename traits::time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> error;
        error.resize(num_nodes + 1);
        std::generate(error.begin(), error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const typename traits::time_t ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), this->get_states()[m]);
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled> //vector encap_ raumdaten an jedem Quadraturpunkt
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                            const typename SweeperTrait::time_t& t)
      {
        UNUSED(t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> rel_error;
        rel_error.resize(error.size());
        std::generate(rel_error.begin(), rel_error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          rel_error[m]->scaled_add(1.0 / this->get_states()[m]->norm0(), error[m]);
        }

        return rel_error;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(),  "evaluating EXPLICIT part at t=" << t);

        auto result = this->get_encap_factory().create();
	auto Mresult = this->get_encap_factory().create();
        result->zero();

	for (int i=0; i<result->data().size(); ++i)
	{result->data()[i]= 8*this->_nu*this->_nu*u->data()[i]*u->data()[i]*(1.00-u->data()[i])/(this->_delta*this->_delta); /*this->source(t)->data()[i]*/;
	//std::cout << u->data()[i] << " "<< result->data()[i] <<" "<< this->_delta <<std::endl; 
	}
	
	this->M_dune.mv(result->data(), Mresult->data());
	
        this->_num_expl_f_evals++;
	
        return Mresult;
	//return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
	
        ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);



        auto result = this->get_encap_factory().create();
        auto rhs = this->get_encap_factory().create();
        auto rhs2 = this->get_encap_factory().create();
        double nu =this->_nu;

	
	this->A_dune.mmv(u->get_data(), result->data());
        /*this->A_dune.mmv(u->get_data(), rhs->data());
        this->A_dune.mmv(u->get_data(), rhs2->data());




        //auto DirichletValues = [] (auto x) {return (x[0]<1e-8 or x[0]>0.9999) ? 0 : x[0];};

        Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(this->M_dune);

        Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(this->M_dune,1.0);

        Dune::CGSolver<VectorType> cg(linearOperator,
                                preconditioner,
                                1e-10, // desired residual reduction factor
                                5000,    // maximum number of iterations
                                0);    // verbosity of the solver
        Dune::InverseOperatorResult statistics ;

        cg.apply(result->data(), rhs->data() , statistics );

        auto var = this->get_encap_factory().create();
        this->M_dune.mv(result->data(), var->data());
        auto neu = this->get_encap_factory().create();
        rhs2->scaled_add(-1.0 , var)  ;*/

        /*result->data() *= nu;
        std::cout << "f_impl mit evaluate" << std::endl;
        for (size_t i = 0; i < result->data().size(); i++) {
          std::cout << result->data()[i] << std::endl;
        }*/


        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {
	

        Dune::BlockVector<Dune::FieldVector<double,1> > M_rhs_dune ;
        M_rhs_dune.resize(rhs->get_data().size());
	
	
	M_rhs_dune = rhs->get_data(); 

        //this->M_dune.mv(rhs->data(), M_rhs_dune); //multipliziert rhs mit matrix_m_dune

	
	
	auto isLeftDirichlet = [] (auto x) {return (x[0] < -20.0 + 1e-8 ) ;};
	auto isRightDirichlet = [] (auto x) {return (x[0] > 20.0 - 1e-8 ) ;};
 
	
	
	std::vector<double> dirichletLeftNodes;
	interpolate(*basis, dirichletLeftNodes, isLeftDirichlet);  

  	std::vector<double> dirichletRightNodes;
  	interpolate(*basis, dirichletRightNodes, isRightDirichlet);
	
  for(int i=0; i<rhs->data().size(); ++i){
  
    if(dirichletLeftNodes[i])
    M_rhs_dune[i] = 1;

        if(dirichletRightNodes[i])
          M_rhs_dune[i] = 0;


    
  }
	
	
	
        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dtA_dune = 	Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->A_dune);
        M_dtA_dune *= (dt * this->_nu);
        M_dtA_dune += this->M_dune;

	auto isDirichlet = [] (auto x) {return (x[0] < -20.0 + 1e-8 or x[0] > 20.0-1e-8) ;};
        std::vector<char> dirichletNodes;
        interpolate(*basis, dirichletNodes, isDirichlet);  //tu valjda interpoliramo kao na nrpdju

        for (size_t i=0; i< M_dtA_dune.N(); i++){
            if (dirichletNodes[i]){
                auto cIt = M_dtA_dune[i].begin();
                auto cEndIt = M_dtA_dune[i].end();
                for(; cIt!=cEndIt; ++cIt){
                    *cIt = (i==cIt.index()) ? 1.0 : 0.0;// 0.0;
                }
            }
        }

        using BitVector = Dune::Solvers::DefaultBitVector_t<VectorType>;
        BitVector ignore(u->data().size());
        for(int i=0; i<rhs->data().size(); ++i){
          if (dirichletNodes[i])
            ignore[i]=true;
        }

        //// Transfer setup
        ////
        using TransferOperator = CompressedMultigridTransfer<VectorType>;
        using TransferOperators = std::vector<std::shared_ptr<TransferOperator>>;

        auto gridptr = FinEl->get_grid();
        TransferOperators transfer(gridptr->maxLevel());
        for (size_t i = 0; i < transfer.size(); ++i)
        {
          // create transfer operator from level i to i+1
          transfer[i] = std::make_shared<TransferOperator>();
          transfer[i]->setup(*gridptr, i, i+1);
        }

        auto& uu = u->data();
        std::cout << "u0 =" << uu[0] << " u[end]=" << uu[uu.size()-1] << std::endl;
       

        //// TNNMG without actual constraints
        VectorType lower(u->data().size());
        lower = -989999;
        VectorType upper(u->data().size());
        upper = 342434;
        using Functional = Dune::TNNMG::BoxConstrainedQuadraticFunctional<MatrixType&, VectorType&, VectorType&, VectorType&, double>;

        auto J = Functional(M_dtA_dune, M_rhs_dune, lower, upper);

        auto localSolver = gaussSeidelLocalSolver(Dune::TNNMG::ScalarObstacleSolver());
        //auto localSolver = gaussSeidelLocalSolver(TrivialSolver());

        using NonlinearSmoother = Dune::TNNMG::NonlinearGSStep<Functional, decltype(localSolver), BitVector>;
        auto nonlinearSmoother = std::make_shared<NonlinearSmoother>(J, u->data(), localSolver);


        using Linearization = Dune::TNNMG::BoxConstrainedQuadraticFunctionalConstrainedLinearization<Functional, BitVector>;
        using DefectProjection = Dune::TNNMG::ObstacleDefectProjection;
        using LineSearchSolver = TrivialSolver;

        using MultiGrid =Dune::Solvers::MultigridStep<MatrixType, VectorType, BitVector>;
        auto mgStep = std::make_shared<MultiGrid>();
        auto gssmoother = Dune::Solvers::BlockGSStepFactory<MatrixType, VectorType, BitVector>::create(Dune::Solvers::BlockGS::LocalSolvers::gs());
        mgStep->setSmoother(&gssmoother);
        mgStep->setTransferOperators(transfer);
        mgStep->setMGType(1,3,3);
        //BitVector dummy(ignore.size());
        //mgStep->setIgnore(dummy);
        //mgStep->setIgnore(ignore);

        auto umfpack = Dune::Solvers::UMFPackSolver<MatrixType, VectorType>{};
        mgStep->basesolver_=&umfpack;

        //std::cout << ENABLE_SUITESPARSE << std::endl;
        //std::cout << HAVE_SUITESPARSE_UMFPACK << std::endl;
        using Step = Dune::TNNMG::TNNMGStep<Functional, BitVector, Linearization, DefectProjection, LineSearchSolver>;
        using Solver = LoopSolver<VectorType>;
        using Norm =  EnergyNorm<MatrixType, VectorType>;

        int mu=1;
        auto projection = DefectProjection();
        auto step = Step(J, u->data(), nonlinearSmoother, mgStep, mu, projection, LineSearchSolver());
        step.setIgnore(ignore);
        //mgStep->setProblem(M_dtA_dune, u->data(), M_rhs_dune);
        auto norm = Norm(M_dtA_dune);
        //auto solver = Solver(&(*mgStep), 1e9, 1e-8, &norm, Solver::FULL);
        auto solver = Solver(&step, 1e9, 0, &norm, Solver::FULL);


        solver.addCriterion(
            [&](){
            return Dune::formatString("   % 12.5e", J(u->data()));
            },
            "   energy      ");

        double initialEnergy = J(u->data());
        solver.addCriterion(
            [&](){
            static double oldEnergy=initialEnergy;
            double currentEnergy = J(u->data());
            double decrease = currentEnergy - oldEnergy;
            oldEnergy = currentEnergy;
            return Dune::formatString("   % 12.5e", decrease);
            },
            "   decrease    ");

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

        std::vector<double> correctionNorms;
        auto tolerance = 1e-8;
        solver.addCriterion(Dune::Solvers::correctionNormCriterion(step, norm, tolerance, correctionNorms));


        solver.preprocess();
        solver.solve();

	
        //Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dtA_dune);

        //Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dtA_dune,1.0);

        //Dune::CGSolver<VectorType> cg(linearOperator,
                              //preconditioner,
                              //1e-10, // desired residual reduction factor
                              //5000,    // maximum number of iterations
                              //0);    // verbosity of the solver

        //Dune::InverseOperatorResult statistics ;

        //cg.apply(u->data(), M_rhs_dune , statistics ); //rhs ist nicht constant!!!!!!!!!



	
	//f = evaluate_rhs_impl(0, u);
	
        ML_CVLOG(4, this->get_logger_id(),
                 "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);


	
	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);
	

        std::cout << "f_impl mit impl_solve" << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          //f->data()[i] = (u->data()[i] - rhs->data()[i]) / (dt);
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
	  //std::cout << "F U " << f->data()[i] << std::endl;
        }
        
        /*for(int i=0; i<f->data().size(); ++i){
	
	  if(dirichletLeftNodes[i])
	  f->data()[i] = 0;

    	  if(dirichletRightNodes[i])
          f->data()[i] = 0;


	  
	}*/
        
        
        //evaluate_rhs_impl(0, u);
        //std::exit(0);

        this->_num_impl_solves++;
        //if (this->_num_impl_solves==5) std::exit(0);


      }
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
