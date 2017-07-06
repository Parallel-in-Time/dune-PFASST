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

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>

#include <iostream>
//#include <c++/4.8/memory>
// dune solvers
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

struct TrivialSolver {
  template<class Vector, class Functional, class BitVector>
    constexpr void operator()(Vector& x, const Functional& f, const BitVector& ignore) const
    { x=1.0;}

};
template<class Vector, class Functional, class BitVector>
struct TrivialLocalSolver {
    constexpr void operator()(Vector& x, const Functional& f, const BitVector& ignore) const
    { x=0.0;}
};



double pi = 3.14159265359;
double two_pi = 2*pi;
double pi_sqr= pi*pi;

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
      Heat_FE<SweeperTrait, Enabled>::Heat_FE(std::shared_ptr<Dune::Functions::PQkNodalBasis<GridType::LevelGridView,SweeperTrait::BASE_ORDER>> basis, size_t nlevel, std::shared_ptr<GridType> grid)
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
	//this->FinEl = FinEl;
	//basis = FinEl->get_basis(nlevel);
	this->nlevel = nlevel;    
	    
	this->basis = basis;

        this->grid = grid;
	
        assembleProblem(basis, this->A_dune, this->M_dune);

        stiffnessMatrix = this->A_dune;
        stiffnessMatrix *= -1;
        w = std::make_shared<VectorType>(this->M_dune.M());
        for(int j=0; j<this->M_dune.M(); j++){
            (*w)[j]=0;
        }

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

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::set_options()
      {


        std::cout << " set_options " <<  std::endl;  
        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        int num_nodes = this->get_quadrature()->get_num_nodes();

        //assembleProblem(basis, A_dune, M_dune);
        std::cout << " set_options " <<  std::endl; 

      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        /*auto result = this->get_encap_factory().create();
        const auto dim = 1; //SweeperTrait::DIM;
	spatial_t n  = this-> _n;
        spatial_t l0 = this-> _nu;
	spatial_t l1 = l0/2. *(pow((1+n/2.), 1/2.) + pow((1+ n/2.), -1/2.) );
	spatial_t d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
        auto exact_solution = [l0, l1, n, d, t](const Dune::FieldVector<double,dim>&x){ 
	  return pow((1 + (pow(2, n/2.)-1 )* exp(-(n/2.)*d*(x+2*l1*t)) ), -2./n);
        };
        auto N_x = [t](const Dune::FieldVector<double,dim>&x){
            return x;
        };
        Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
        interpolate(*basis, x_node, N_x);
        interpolate(*basis, result->data(), exact_solution);
        return result;*/
        
        auto result = this->get_encap_factory().create();

        const auto dim = 2;//SweeperTrait::DIM;
        spatial_t nu = this-> _nu; 
	
	
	auto exact_solution1 = [t,  nu, dim](const Dune::FieldVector<double,dim>&x){
             
	   
	   double eps = 1e-8;
	   
	  if( (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) < 0.25)
	    return 0.5;
	    return 1.0;
            
        };

	auto exact_solution2 = [t, nu, dim](const Dune::FieldVector<double,dim>&x){
          double eps = 1e-8; 

	   if( (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) < 0.0025)
	   return 0.25;
	   return 0.0;
	};
	
	 auto N_x = [t](const Dune::FieldVector<double,dim>&x){
            return x;

        };

        Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
	Dune::BlockVector<Dune::FieldVector<double,1>> pom1, pom2;
	pom1.resize(basis->size());
	pom2.resize(basis->size());
        
	interpolate(*basis, x_node, N_x);

	interpolate(*basis, pom1, exact_solution1);
	interpolate(*basis, pom2, exact_solution2);
	
	
	
	for (int i = 0; i< basis->size(); ++i)
	{
	  result->data()[i][0] = pom1[i];
	  //result->data()[i][1] = pom2[i];
	  
	  
	}

        //auto grid = (*(this->get_fine())).get_grid();
        typedef Dune::BlockVector<Dune::FieldVector<double, 2> > VectorType;
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ColumnType;
        typedef Dune::YaspGrid<2> GridType; //ruth_dim
        //typedef GridType::LevelGridView GridView;
        //GridType::LevelGridView gridView = grid->leafGridView();
        //Dune::VTKWriter<GridView> vtkWriter(gridView);



        /*string name = "initial";
        ColumnType sol_u, sol_v;
        sol_u.resize(result->data().size());
        sol_v.resize(result->data().size());
        for (int i =0; i< result->data().size(); ++i)
        {
          sol_u[i] = result->data()[i][0];
          sol_v[i] = result->data()[i][1];
        }
        vtkWriter.addVertexData(sol_u, "fe_solution_u");
        vtkWriter.addVertexData(sol_v, "fe_solution_v");
        vtkWriter.write("gray_scott" + name);*/




return result;
        
        
      }
      
      
      
      
      //______________________________
      
      
      
      /*template<class SweeperTrait, typename Enabled>
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
      }*/
      
      
      
      
      
      
      
      
      
      //_______________________________
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        //ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        //this->_num_expl_f_evals = 0;
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

      //typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; //ruth_dim
     
     
      
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

      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(),  "evaluating EXPLICIT part at t=" << t);

        auto result2 = this->get_encap_factory().create();
	auto Mresult = this->get_encap_factory().create();
        result2->zero();

	
	
        this->_num_expl_f_evals++;
	
        return result2;

      }*/

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
	
         ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);



        auto result = this->get_encap_factory().create();
        auto u2 = this->get_encap_factory().create();
        double nu =this->_nu;

	u2->zero();
	for (int i=0; i<u->get_data().size(); ++i)
	    {
	    result->data()[i]= -pow(u->get_data()[i], _n+1) * (*w)[i];	
	    }
	//this->M_dune.mv(u2->get_data(), result->data());
	this->M_dune.umv(u->get_data(), result->data());
	result->data()*=_nu*_nu;
	this->A_dune.umv(u->get_data(), result->data());

	

        //result->data() *= nu;
	/*std::cout << "evaluate  " << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {

	  std::cout << "f " << result->data()[i] << std::endl;

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
	


 int my_rank;  
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );  
        //MPI_Barrier(MPI_COMM_WORLD);

        std::cout << my_rank << " ***ANFANG u: " << u->norm0() << "num_level" << this->is_coarse << std::endl;
	
        ML_CVLOG(4, this->get_logger_id(),
                 "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);



	

    // solve implicit through TNNMG (actually: Newton)
    {

	  auto isLeftDirichlet = [] (auto x) {return (x[0] < -200 + 1e-8 ) ;};
	  auto isRightDirichlet = [] (auto x) {return (x[0] > 200 - 1e-8 ) ;};
 
	
	
	  std::vector<double> dirichletLeftNodes;
	  interpolate(*basis, dirichletLeftNodes, isLeftDirichlet);  

	  std::vector<double> dirichletRightNodes;
	  interpolate(*basis, dirichletRightNodes, isRightDirichlet);


          auto phiprime =
            [&](auto&& uu) {
            return dt*_nu*_nu*pow(uu, _n+1);
          };

          auto phi2prime = 
            [&] (auto&& uu) {
            return dt*(_nu*_nu)*(_n+1) * pow(uu, _n);
          };

	
          auto matrix_ = this->M_dune;
          matrix_ *= (1-dt*_nu*_nu);
          matrix_.axpy(-dt, this->A_dune);
  
          //std::cout << "Anfang TNNMG" << std::endl;
          using BitVector = Dune::Solvers::DefaultBitVector_t<VectorType>;
          BitVector ignore(u->data().size());

          //// Transfer setup
          using TransferOperator = CompressedMultigridTransfer<VectorType>;
          using TransferOperators = std::vector<std::shared_ptr<TransferOperator>>;


        auto& gridptr = this->grid;
        int num_level=0;
        if(nlevel==0 ){
            num_level= gridptr->maxLevel();
            
        }else{
            num_level= gridptr->maxLevel() -1; 
        }
        TransferOperators transfer(num_level);
        for (size_t i = 0; i < transfer.size(); ++i)
        {
          // create transfer operator from level i to i+1
          transfer[i] = std::make_shared<TransferOperator>();
          transfer[i]->setup(*gridptr, i, i+1);
        }
        
        //std::cout << "impl solve TNNMG u: " << u->get_data()[u->get_data().size()-1] << "num_level" << this->is_coarse << std::endl;

        auto& uu = u->data();
        //std::cout << "u0 =" << uu[0] << " u[end]=" << uu[uu.size()-1] << "numlevel " << num_level << std::endl;


        //// TNNMG without actual constraints
        //using Functional = Dune::TNNMG::EnergyFunctional<MatrixType, VectorType, double,1>;
        using Functional = Dune::TNNMG::EnergyFunctional<MatrixType, VectorType, decltype(phiprime), decltype(phiprime), decltype(phi2prime), double>;

        //auto J = Functional(df, newton_rhs, lower, upper);
        //auto J = Functional(this->M_dune, this->A_dune, dt, _nu, *(this->w), rhs->data());
        auto J = Functional(matrix_, rhs->data(), *(this->w), phiprime, phiprime, phi2prime);
        //using LocalFunctional = Dune::TNNMG::EnergyFunctional<MatrixType::block_type, VectorType::block_type, decltype(phiprime), decltype(phiprime), decltype(phi2prime), double>;
        //using LocalFunctional = Dune::TNNMG::EnergyDirectionalRestriction<MatrixType::block_type, VectorType::block_type, decltype(phiprime),double>;


        //auto localSolver = gaussSeidelLocalSolver(Dune::TNNMG::ScalarObstacleSolver());
        auto localSolver = gaussSeidelLocalSolver(Dune::TNNMG::ScalarBisectionSolver());
        //auto localSolver = Dune::TNNMG::ScalarBisectionSolver();   //
 
        //auto localSolver = gaussSeidelLocalSolver(TrivialLocalSolver<double,Dune::TNNMG::EnergyFunctional<double,double, DUMMYPHI, decltype(phiprime), decltype(phi2prime), double>, BitVector::value_type>()); // trivial solver does not change anything
        //auto localSolver = TrivialLocalSolver<double,LocalFunctional, BitVector::value_type>(); // trivial solver does not change anything
        //auto localSolver = gaussSeidelLocalSolver(TrivialLocalSolver<double,LocalFunctional, BitVector::value_type>()); // trivial solver does not change anything

        using NonlinearSmoother = Dune::TNNMG::NonlinearGSStep<Functional, decltype(localSolver), BitVector>;
        auto nonlinearSmoother = std::make_shared<NonlinearSmoother>(J, u->data(), localSolver);


        //using Linearization = Dune::TNNMG::BoxConstrainedQuadraticFunctionalConstrainedLinearization<Functional, BitVector>;
        using Linearization = Dune::TNNMG::EnergyFunctionalConstrainedLinearization<Functional, BitVector>; // incorporates what used to be evaluate_f, evaluate_df
        //using DefectProjection = Dune::TNNMG::ObstacleDefectProjection;
        //using LineSearchSolver = TrivialSolver; // always gives 1 as correction damping
        using LineSearchSolver = Dune::TNNMG::ScalarBisectionSolver; // always gives 1 as correction damping
        //using LineSearchSolver = Dune::TNNMG::ScalarObstacleSolver;

        /* Setup linear multigrid */
        using MultiGrid =Dune::Solvers::MultigridStep<MatrixType, VectorType, BitVector>;
        auto mgStep = std::make_shared<MultiGrid>();
        auto gssmoother = Dune::Solvers::BlockGSStepFactory<MatrixType, VectorType, BitVector>::create(Dune::Solvers::BlockGS::LocalSolvers::gs());
        mgStep->setSmoother(&gssmoother);
        mgStep->setTransferOperators(transfer);
        mgStep->setMGType(1,3,3);

        // base solver for multigrid
        auto umfpack = Dune::Solvers::UMFPackSolver<MatrixType, VectorType>{}; // direct solver
        mgStep->basesolver_=&umfpack;


        auto trivialProjection = [](auto&& f, auto& x, auto& c) {};
        //using Step = Dune::TNNMG::TNNMGStep<Functional, BitVector, Linearization, DefectProjection, LineSearchSolver>;
        using Step = Dune::TNNMG::TNNMGStep<Functional, BitVector, Linearization, decltype(trivialProjection), LineSearchSolver>;
        int mu=1;
        //
        // J Functional to minimize, mgStep: linear "solver" for linear correction, mu #multigrid steps per iteration, trivialProjection is identity
        auto step = Step(J, u->data(), nonlinearSmoother, mgStep, mu, trivialProjection, LineSearchSolver());
        //step.setPreSmoothingSteps(0); nach lasse reintun
        step.setIgnore(ignore);

        //using Norm =  TwoNorm<VectorType>;
        using Norm =  EnergyNorm<MatrixType,VectorType>;
        auto norm = Norm(stiffnessMatrix);
        // max. 20 iterations or two norm of correction less than 1e-10
        using Solver = LoopSolver<VectorType>;
        auto solver = Solver(&step, 20, 1e-10, &norm, Solver::FULL);


        //solver.addCriterion(
            //[&](){
            //return Dune::formatString("   % 12.5e", J(u->data()));
            //},
            //"   energy      ");

        //double initialEnergy = J(u->data());
        //solver.addCriterion(
            //[&](){
            //static double oldEnergy=initialEnergy;
            //double currentEnergy = J(u->data());
            //double decrease = currentEnergy - oldEnergy;
            //oldEnergy = currentEnergy;
            //return Dune::formatString("   % 12.5e", decrease);
            //},
            //"   decrease    ");

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

        //std::vector<double> correctionNorms;
        //auto tolerance = 1e-8;
        //solver.addCriterion(Dune::Solvers::correctionNormCriterion(step, norm, tolerance, correctionNorms));

        solver.preprocess();

        solver.solve();


  
	}

	
	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);

    //std::cout << "impl solve "  << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
	  //std::cout << "f " << f->data()[i] << std::endl;
        }
        //evaluate_rhs_impl(0, u);
	//std::exit(0);
        this->_num_impl_solves++;
        //if (this->_num_impl_solves==5) std::exit(0);
        std::cout << my_rank << " ENDE u: " << u->norm0() << " " <<  this->is_coarse<< std::endl;
        
        //std::exit(0);




      }
      
      template<class SweeperTrait, typename Enabled>
        void
        Heat_FE<SweeperTrait, Enabled>::evaluate_f(shared_ptr<typename SweeperTrait::encap_t> f,
            const shared_ptr<typename SweeperTrait::encap_t> u,
            const typename SweeperTrait::time_t& dt,
            const shared_ptr<typename SweeperTrait::encap_t> rhs
            ){


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

	
	/*f->zero();
	
	Dune::BlockVector<Dune::FieldVector<double,1> > fneu;
        fneu.resize(u->get_data().size());
	for (int i=0; i<f->data().size(); ++i)
	{fneu[i]= -8*this->_nu*this->_nu *u->data()[i]*u->data()[i]*(1.00-u->data()[i])/(this->_delta*this->_delta);

	}
	
	this->M_dune.mv(fneu, f->data());
	
	this->A_dune.mmv(u->get_data(),f->data());
	f->data() *= dt;
	this->M_dune.umv(u->get_data(),f->data());
	f->data() -=rhs->get_data();*/

      }
						
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						){
          
          
          
          
          
            for (int i=0; i<df.N(); ++i)
            {
                df[i][i]= (_nu*_nu)*(_n+1) * pow(u->get_data()[i], _n) * ((double) (*w)[i]);	
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
          
          
          


      }  
      
      
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
