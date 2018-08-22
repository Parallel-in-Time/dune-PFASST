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




#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/densematrix.hh>
#include<dune/common/parallel/mpihelper.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/io.hh>
#include<dune/istl/matrixmarket.hh>
#include<dune/istl/matrixredistribute.hh>
#include <dune/istl/schwarz.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/typetree/utility.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#if USE_DG
#  include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#else
#  include <dune/functions/functionspacebases/pqknodalbasis.hh>
#endif
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/istlbackend.hh>
#include <dune/fufem/formatstring.hh>
#include <dune/fufem/assemblers/transferoperatorassembler.hh>

#include<dune/istl/paamg/pinfo.hh>
#include<dune/istl/paamg/graph.hh>
#if USE_DG
#  include <dune/parmg/test/dglaplacematrix.hh>
#else
#  include <dune/parmg/test/laplacematrix.hh>
#endif
#include <dune/parmg/iterationstep/lambdastep.hh>
#include <dune/parmg/iterationstep/multigrid.hh>
#include <dune/parmg/iterationstep/multigridstep.hh>
#include <dune/parmg/norms/normadapter.hh>
#include <dune/parmg/parallel/communicationp1.hh>
#include <dune/parmg/parallel/communicationdg.hh>
#include <dune/parmg/parallel/datahandle.hh>
#include <dune/parmg/parallel/dofmap.hh>
#include <dune/parmg/parallel/globaldofindex.hh>
#include <dune/parmg/parallel/istlcommunication.hh>
#include <dune/parmg/parallel/matrixalgebra.hh>
#include <dune/parmg/parallel/vectoralgebra.hh>
#include <dune/parmg/parallel/redistributematrix.hh>
#include <dune/parmg/parallel/redistributevector.hh>
#include <dune/parmg/parallel/parallelenergyfunctional.hh>
#include <dune/parmg/parallel/parallelenergynorm.hh>
#include <dune/parmg/parallel/restrictmatrix.hh>
#include <dune/parmg/solvers/coarsesuperlusolver.hh>
#include <dune/parmg/solvers/linesearch.hh>
#include <dune/parmg/solvers/directionsearch.hh>
#include <dune/parmg/iterationstep/multigridsetup.hh>
#include <dune/parmg/iterationstep/parallelprojectedgs.hh>

#include <dune/solvers/iterationsteps/blockgssteps.hh>
#include <dune/solvers/norms/energynorm.hh>
#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/transferoperators/compressedmultigridtransfer.hh>



using namespace std;
//using namespace Dune;

//using namespace Dune;
using namespace Dune::ParMG;

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
	std::cout << "get basis " << std::endl;
	basis = FinEl->get_basis(nlevel);
	this->nlevel=nlevel;
	grid = FinEl->get_grid();
	    
	//assembleProblem(basis, this->A_dune, this->M_dune);
	auto sBackend = Dune::Fufem::istlMatrixBackend(this->A_dune);
        auto mBackend = Dune::Fufem::istlMatrixBackend(this->M_dune);
	
        using Assembler = Dune::Fufem::DuneFunctionsOperatorAssembler<BasisFunction, BasisFunction>;
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

        const auto bs = basis->size();
        std::cout << "Finite Element basis of level " << nlevel << " consists of " <<  basis->size() << " elements " << std::endl;

	int rank, num_pro;
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank );
    	MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
    	

    	
  	/*using MGSetup = Dune::ParMG::ParallelMultiGridSetup< BasisFunction, MatrixType, VectorType >;
  	MGSetup mgSetup{*grid,grid->maxLevel() - (this->nlevel)};
  	auto gridView = mgSetup.bases_.back().gridView();
 	using MG = Dune::ParMG::Multigrid<VectorType>;
  	MG mg;
  
    	using namespace Dune::ParMG;
    	auto& levelOp = mgSetup.levelOps_;
    	auto df_pointer = std::make_shared<MatrixType>(this->A_dune);	
    	mgSetup.matrix(df_pointer);
    	auto fineIgnore = std::make_shared< Dune::BitSetVector<1> >(bs);
    	for (std::size_t i = 0; i < bs; ++i){
      		(*fineIgnore)[i] = false;
      		if(i==0 &&rank==0) (*fineIgnore)[i] = true;
      		if(rank==num_pro - 1 && i== bs-1) (*fineIgnore)[i] = true;
      	}
      	    		std::cout << "nach ignore gesetzt " << std::endl;

    	mgSetup.ignore(fineIgnore);
    	mgSetup.setupLevelOps();
    	double dampening =1.0;
    	mgSetup.setupSmoother(dampening);
    	bool enableCoarseCorrection=true;
    	if (enableCoarseCorrection)
      		mgSetup.setupCoarseSuperLUSolver();
    	else
      		mgSetup.setupCoarseNullSolver();
    	mg.levelOperations(levelOp);
    	mg.coarseSolver(mgSetup.coarseSolver());
    	//levelOp.back().maybeRestrictToMaster(newton_rhs);
    	std::function<void(VectorType&)> collect = Dune::ParMG::makeCollect<VectorType>(*mgSetup.comms_.back());
    	std::function<void(VectorType&)> restrictToMaster = [op=levelOp.back()](VectorType& x) { op.maybeRestrictToMaster(x); };
    	std::cout << "im sweeper vor energyfunctional vor" << std::endl;
    	

      	auto vz =  this->A_dune;
	vz *=-1;    

		
        auto  parallel_energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(vz, restrictToMaster, gridView.grid().comm());*/
        //std::cout << typeid(parallel_energyNorm).name() << std::endl; std::exit(0);
	//int i= parallel_energyNorm;

        /*int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );

        	MPI_Comm comm_x, comm_t; 
	int myid, xcolor, tcolor;

	int space_num=2;
   	xcolor = (my_rank / space_num);
   	tcolor = my_rank % space_num;

   	MPI_Comm_split( MPI_COMM_WORLD, xcolor, my_rank, &comm_x );
   	MPI_Comm_split( MPI_COMM_WORLD, tcolor, my_rank, &comm_t );*/


        this->encap_factory()->set_size(bs);
        //this->encap_factory()->set_comm(comm_x);
        //this->encap_factory()->set_FE_manager(FinEl, (this->nlevel), this->A_dune);
        this->A_dune*=-1;
        //std::cout << "############################################### groesse mass matrix " << this->A_dune.N() << std::endl;
        //std::cout << "alles assembliert" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);

      }
      

       
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::set_options()
      {


        //std::cout << " set_options " <<  std::endl;  
        IMEX<SweeperTrait, Enabled>::set_options();

        //this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        int num_nodes = this->get_quadrature()->get_num_nodes();

        //assembleProblem(basis, A_dune, M_dune);
        //std::cout << " set_options " <<  std::endl; 
        //MPI_Barrier(MPI_COMM_WORLD);

      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        
        const auto dim = 1; 
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
    	//int rank, num_pro;
    	//MPI_Comm_rank(MPI_COMM_WORLD, &rank );
        interpolate(*basis, result->data(), exact_solution);
        /*for(int i=0; i< result->get_data().size(); i++){
           if(rank==0) std::cout<< rank << " initial " << i << " " << result->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);
        for(int i=0; i< result->get_data().size(); i++){
           if(rank==1) std::cout<< rank << " initial " << " " << i << " " << result->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);*/
        
        
        /*int rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
	

        
	for (size_t i = 0; i < result->get_data().size(); i++) {
	  if (rank==0) std::cout << i << t << " u " << result->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);*/	//std::exit(0);
        

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


      /*template<class SweeperTrait, typename Enabled> //Fehler der aktuellen Loesung an jedem Quadraturpunkt
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
      }*/

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

	//std::cout << "++++++++++++ " << this->_nu << std::endl;
	u2->zero();
	for (int i=0; i<u->get_data().size(); ++i)
	    {
	    u2->data()[i]= -pow(u->get_data()[i], _n+1);	
	    }
	this->M_dune.mv(u2->get_data(), result->data());
	this->M_dune.umv(u->get_data(), result->data());
	result->data()*=(this->_nu*this->_nu);
	this->A_dune.umv(u->get_data(), result->data());

	//std::cout << "im evaluate " << std::endl;
	//std::exit(0);

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



////////////////////////////////////////////////////////////////////////////// non linear cg                 
       /*              auto residuum = this->get_encap_factory().create();
	Dune::BlockVector<Dune::FieldVector<double,1> > newton_rhs, newton_rhs2 ;
    newton_rhs.resize(rhs->get_data().size());
    newton_rhs2.resize(rhs->get_data().size());



	for (int i=0; i< 200 ;i++){
	  Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune); ///////M
	  evaluate_f(f, u, dt, rhs);
	  evaluate_df(df, u, dt);
	  df.mv(u->data(), newton_rhs);
	  newton_rhs -= f->data();
      	  newton_rhs2 = newton_rhs; 
	

  
	  Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df);
	  Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df,1);
	  Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-12, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver
	  Dune::InverseOperatorResult statistics ;
	  cg.apply(u->data(), newton_rhs , statistics );*/


////////////////////////////////////////////////////////////////////////////////////////// ende sequential start parallel solever

	/*for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==0) std::cout << i << "_rhs " << rhs->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
        for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==num_pro-1) std::cout << i << " rhs " << rhs->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);std::exit(0);*/

	int rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
	

        auto residuum = this->get_encap_factory().create();
	Dune::BlockVector<Dune::FieldVector<double,1> > newton_rhs, newton_rhs2 ;
        newton_rhs.resize(rhs->get_data().size());
        newton_rhs2.resize(rhs->get_data().size());
        

        auto delta_u = this->get_encap_factory().create(); //u->data();
        delta_u->data()*=0;

        
	for (int i=0; i< 10 ;i++){
	  Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune); ///////M

	  
	  
	  evaluate_f(f, u, dt, rhs);

	  evaluate_df(df, u, dt);

	  newton_rhs*=0;
	  newton_rhs -= f->data();
      	  newton_rhs2 = newton_rhs;
      	  

	//std::cout  << rank << " im impl solve" << std::endl;
     
	
	///////////////////////////////////////////////////////////////////////////////////////

	//here we use parmg solver to solve the local system (df * u = newton_rhs) for the unknown u
  	//auto &x=delta_u->data();//abst
  	using MGSetup = Dune::ParMG::ParallelMultiGridSetup< BasisFunction, MatrixType, VectorType >;


  	MGSetup mgSetup{*grid, grid->maxLevel() - (this->nlevel)};//sdc -1 mlsdc - (this->nlevel)}; //0 grobgitter

  	auto gridView = mgSetup.bases_.back().gridView();
 	using MG = Dune::ParMG::Multigrid<VectorType>;
  	MG mg;
        //if (this->get_communicator()->get_rank()!=rank) std::exit(0);
    	using namespace Dune::ParMG;
    	auto& levelOp = mgSetup.levelOps_;
    	
    	/*auto isDirichlet_left = [] (auto x) {return (x[0] < -200.0 + 1e-8 ) ;};
    	auto isDirichlet_right = [] (auto x) {return (x[0] > 200.0-1e-8) ;};
  	std::vector<char> dirichletNodes_right, dirichletNodes_left;
  	interpolate(*basis, dirichletNodes_right, isDirichlet_right); 
  	interpolate(*basis, dirichletNodes_left, isDirichlet_left); */
  	
  	/*for (size_t i=0; i<df.N(); i++){
    		if (dirichletNodes_left[i]||dirichletNodes_right[i]){
      			auto cIt = df[i].begin();
      			auto cEndIt = df[i].end();
      		for(; cIt!=cEndIt; ++cIt){
        		*cIt = (i==cIt.index()) ? 1.0 : 0.0; // 0.0;
      		}
      		if(dirichletNodes_left[i]) newton_rhs[i]=0;
      		else newton_rhs[i]=0;
    		}
  	}*/
    	
    	/*MPI_Barrier(MPI_COMM_WORLD);
    	  for (int i=0; i<df.N(); ++i)
            {
            for (int j=0; j<df.M(); ++j)
                {
                    if (df.exists(i,j)) {
                        if(rank==0) std::cout << df[i][j] << " " ;}	
                }
            }std::cout<<std::endl;*/
 
    	
    	auto df_pointer = std::make_shared<MatrixType>(df);	
    	mgSetup.matrix(df_pointer);
    	auto fineIgnore = std::make_shared< Dune::BitSetVector<1> >(u->get_data().size());
    	for (std::size_t i = 0; i < u->get_data().size(); ++i){
      		(*fineIgnore)[i] = false;
      		//if(i==0&&(rank==0) ) (*fineIgnore)[i] = true; //
      		//if((rank==num_pro - 1)&&i== u->get_data().size()-1) (*fineIgnore)[i] = true; //
      	}
      	    		//std::cout << "nach ignore gesetzt " << std::endl;
      	//std::cout << "im sweeper ignore gestzt" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
    	mgSetup.ignore(fineIgnore);
    	mgSetup.setupLevelOps();
    	double dampening =1.0;
    	mgSetup.setupSmoother(dampening);
    	bool enableCoarseCorrection=true;
    	if (enableCoarseCorrection)
      		mgSetup.setupCoarseSuperLUSolver();
    	else
      		mgSetup.setupCoarseNullSolver();
    	mg.levelOperations(levelOp);
    	mg.coarseSolver(mgSetup.coarseSolver());
    	levelOp.back().maybeRestrictToMaster(newton_rhs);
    	std::function<void(VectorType&)> collect = Dune::ParMG::makeCollect<VectorType>(*mgSetup.comms_.back());
    	std::function<void(VectorType&)> restrictToMaster = [op=levelOp.back()](VectorType& x) { op.maybeRestrictToMaster(x); };
    	//std::cout << "im sweeper vor energyfunctional" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	
    	auto energyFunctional = Dune::ParMG::makeParallelEnergyFunctional(
      		*df_pointer,
      		newton_rhs,
      		gridView.grid().comm(),
      		//collect
      		restrictToMaster
      	);
      	
      	//auto vz =  this->A_dune;
	//vz *=-1;
    	auto energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(*df_pointer, restrictToMaster, gridView.grid().comm()); //*df_pointer
    	levelOp.back().maybeCopyFromMaster(delta_u->data()); //abst
    	double tol = 1e-8;
    	//std::cout << "im sweeper vor apply" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	
        collect(newton_rhs);
    	auto realIterationStep = [&](VectorType& x) {
    	        //std::cout << "******************** realIterationStep " << std::endl;
      		auto b = newton_rhs;
      		//std::cout << "******************** realIterationStep 1" << std::endl;
      		mg.apply(x, b);
      		//std::cout << "******************** realIterationStep 2" << std::endl;
    	};
  	auto& feBasis = mgSetup.bases_.back();
    	auto solverNorm = std::make_shared< NormAdapter<VectorType> >(energyNorm);
	auto iterationStep = std::make_shared< LambdaStep<VectorType> >(realIterationStep, delta_u->data()); //abst
	int steps = 500;
    	auto solver = Dune::Solvers::LoopSolver<VectorType>(iterationStep, steps, tol, solverNorm,NumProc::QUIET);// NumProc::FULL); //
    		//std::cout << "im sweeper vorm solver aufruf " << std::endl;
	solver.preprocess();
	
	//solver.setOption(UMFPACK_PRL, 0);
	//std::cout  << rank << " " << this->_num_impl_solves << " vor solver solve" << std::endl;
    	solver.solve();
    	//std::cout  << rank << " " << this->_num_impl_solves << " nach solver solve" << dt << std::endl;
    	u->data()+=delta_u->data();
    	num_solves++;

//////////////////////////////////////////////////////////////////////////////////////
	

	
           //evaluate_f(f, u, rhs, M_dune, A_dune, *w); 


	  	/*for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==0) std::cout << "newton_rhs " << u->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
        for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==1) std::cout << "newton_rhs " << u->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);//std::exit(0);*/
	/*for (int i=0; i< u->data().size(); i++){
          if(rank==0) std::cout << rank << " " << i << " " <<  u->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
        for (int i=0; i< u->data().size(); i++){
          if(rank==1) std::cout  << rank << " " << i << " " << u->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);std::exit(0);*/
           
           	  /*for (size_t i = 0; i < u->get_data().size(); i++) {

	  if (rank ==0) std::cout <<i << "newton  u " << u->data()[i] << std::endl;

        } 	MPI_Barrier(MPI_COMM_WORLD); std::exit(0);*/
        
           evaluate_f(f, u, dt, rhs);            
                     
           ///////////////////////////////////////////////////using Norm =  EnergyNorm<MatrixType,VectorType>;
  
	    
//MPI_Barrier(MPI_COMM_WORLD); std::exit(0);
	   //auto parallel_energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(vz, restrictToMaster, gridView.grid().comm());
	    
            //if (rank==0) std::cout << i << " residuumsnorm von f_global(u) infinity " << norm_global(f) << std::endl;  	                  
            //if (rank==0) std::cout << i << " residuumsnorm von f_global(u) energy " << f.infinity_norm() << std::endl;  
           //std::cout << i << " rank "<< rank << " ################################################################################                          residuumsnorm von f(u) " << parallel_energyNorm(f->data()) << std::endl; 
           if(i == 3) break; 
           //std::cout << i << " rank "<< rank << " ################################################################################                          vor energir norm "  << std::endl; 
           //if( parallel_energyNorm(f->data()) < 1e-10) {           std::cout << i << " process rank breaks inner newton " << rank << std::endl;  break;} 
           
           //if( f->norm0() < 1e-10) {           std::cout << i << " process rank breaks inner newton " << rank << std::endl;  break;} 
           
          // std::cout << i << " rank "<< rank << " ################################################################################                          vor energir norm "  << std::endl;  
        /*for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==0) std::cout << i << " u " << u->data()[i] << std::endl;
        }*///MPI_Barrier(MPI_COMM_WORLD);std::exit(0);
  
	  /*Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df);
	  Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df,1);
	  Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-12, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver
	  Dune::InverseOperatorResult statistics ;
	  cg.apply(u->data(), newton_rhs , statistics ); //rhs ist nicht constant!!!!!!!!!
	  num_solves++;

          evaluate_f(f, u, dt, rhs);
          int my_rank, num_pro;

          if(f->norm0()<this->newton){ if(!this->is_coarse) std::cout << my_rank << "***************************************** anzahl iterationen innerer newton " << i+1 << " " << num_solves <<std::endl;   break;}*/
	  

	}

	
	
	//std::exit(0);
	

	
	
	
	
	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);

        //std::cout << "impl solve "  << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
	  //std::cout << "f " << f->data()[i] << std::endl;
        }
        //evaluate_rhs_impl(0, u);
        //evaluate_f(f, u, dt, rhs);
	//std::cout << "***************************************** das neue f " << f->norm0() << std::endl;
	//std::exit(0);
        this->_num_impl_solves++;
        //if (this->_num_impl_solves==5) std::exit(0);
	//std::cout  << rank << " " << this->_num_impl_solves << " ende impl solve" << std::endl;

      }
      
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::evaluate_f(shared_ptr<typename SweeperTrait::encap_t> f,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt,
						 const shared_ptr<typename SweeperTrait::encap_t> rhs
						){
          
          	int rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank );
          f->zero();
	//std::cout << "-------------- "<< this->_nu << std::endl;
	Dune::BlockVector<Dune::FieldVector<double,1> > fneu;
        fneu.resize(u->get_data().size());
	for (int i=0; i<u->get_data().size(); ++i)
	{
	  fneu[i]= (this->_nu*this->_nu) * (pow(u->get_data()[i], _n+1) - u->get_data()[i]);	
	}
	

	//f->data() *= (this->_nu*this->_nu);
	this->M_dune.mv(fneu, f->data());
	this->A_dune.mmv(u->get_data(),f->data());
	f->data() *= dt;
	this->M_dune.umv(u->get_data(),f->data());
	/*for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==0) std::cout << "fneu " << f->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
        for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==1) std::cout << "fneu " << f->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);//std::exit(0);*/ 
	
	f->data() -=rhs->get_data();

	/*for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==0) std::cout << "fneu " << f->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
        for (size_t i = 0; i < u->get_data().size(); i++) {
	  if (rank==1) std::cout << "fneu " << f->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);//std::exit(0);*/          

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
            double _nu=this->_nu;
            double _n=this->_n;
            for (int i=0; i<df.N(); ++i)
            {
		for(int j=0; j< df.M(); j++)
                	if (df.exists(i,j)) df[i][j]= (_nu*_nu)*(_n+1) * pow(u->get_data()[j], _n) * this->M_dune[i][j];	
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
      } 
						
      /*template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						){
          
          
          //std::cout << "############### " << this->_nu<< std::endl;
          
          
          for (int i=0; i<df.N(); ++i)
            {
            for (int j=0; j<df.M(); ++j)
                {
                    if (df.exists(i,j)) 
                        df[i][j]= -(_n+1) *this->M_dune[i][j] * pow(u->get_data()[j], _n);	
                }
            }
            df += this->M_dune;
          for (int i=0; i<df.N(); ++i)
            {
            for (int j=0; j<df.M(); ++j)
                {
                    if (df.exists(i,j)) 
                        df[i][j]= df[i][j]*(- this->_nu*this->_nu) ;	
                }
            }
	    //df *= (- this->_nu*this->_nu);	
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
          
          
          


      } */ 
      
      
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
