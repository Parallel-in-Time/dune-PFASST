#include <memory>
#include <iostream>
#include <vector>

#include "config.h"


// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/newton/newton.hh>



#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>

#include <pfasst/config.hpp>


// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
//#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/newton/newton.hh>


#include "../../datatypes/pdelab_vec.hpp"

#include <dune/grid/common/mcmgmapper.hh> // mapper class
//#include <dune/grid/utility/gridtype.hh>
#include "FE_sweeper.hpp"

#include "spectral_transfer.hpp"



/*template<typename Number>
class Problem
{
public:
  typedef Number value_type;
  const double dt=0;	
  Problem () {}
  Problem<Number> (double dt_): dt(dt_) {}
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=1;
	const double t=0;	
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
        return solution * std::exp(-t * dim * PI*PI);
  }
};*/


	//setup the grid
        //const int dim=DIM;
	const int degree =1;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        //Dune::FieldVector<double,dim> h = {1};	      
	//std::array<int,dim> n;
	//std::fill(n.begin(), n.end(), nelements);
        //std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(h,n));

        //gridp->refineOptions(false); // keep overlap in cells
        //gridp->globalRefine(1);
        typedef Grid::LeafGridView GV;
        //GV gv=gridp->leafGridView();
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        //FEM fem(gv);

  	// Make grid function space
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	//GFS gfs(gv,fem);


	typedef double RF; 
  	//Problem<RF> mass_problem(0.0);
  	// Assemble constraints
  	typedef typename GFS::template
    	ConstraintsContainer<RF>::Type CC;

  	
	// Make a local operator
  	typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;

  	// Make a global operator
  	typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> MatrixType;



using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, GFS, MatrixType>;

using namespace std;
using namespace pfasst::examples::heat_FE;


namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      using sweeper_t = Heat_FE<sweeper_traits<encap_traits_t>, MatrixType>; //das muss matrixtyp sein 
      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::SDC;
      using pfasst::quadrature::QuadratureType;
      using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

      /*shared_ptr<heat_FE_sdc_t> run_sdc(const size_t nelements, const size_t basisorder, const size_t dim, const size_t nnodes,
                                       const QuadratureType& quad_type, const double& t_0,
                                       const double& dt, const double& t_end, const size_t niter)
      {
        

        return sdc;
      }*/
    }
  }
}




  int main(int argc, char** argv) {
    using pfasst::config::get_value;
    using pfasst::quadrature::QuadratureType;
    using pfasst::examples::heat_FE::Heat_FE;

    using sweeper_t = Heat_FE<pfasst::sweeper_traits<encap_traits_t>, MatrixType>;


	Dune::MPIHelper&
      	helper = Dune::MPIHelper::instance(argc, argv);

    //Dune::MPIHelper::instance(argc, argv);
    pfasst::init(argc, argv, sweeper_t::init_opts);

    auto const   nelements = pfasst::config::get_value<size_t>("num_elements", 100); 
    const size_t nnodes    = get_value<size_t>("num_nodes", 3);
    const QuadratureType quad_type = QuadratureType::GaussRadau;
    const double t_0 = 0.0;
    const double dt = get_value<double>("dt", 0.05);
    double t_end = get_value<double>("tend", 0.1);
    size_t nsteps = get_value<size_t>("num_steps", 0);
    if (t_end == -1 && nsteps == 0) {
      ML_CLOG(ERROR, "USER", "Either t_end or num_steps must be specified.");
      throw std::runtime_error("either t_end or num_steps must be specified");
    } else if (t_end != -1 && nsteps != 0) {
      if (!pfasst::almost_equal(t_0 + nsteps * dt, t_end)) {
        ML_CLOG(ERROR, "USER", "t_0 + nsteps * dt != t_end ("
                               << t_0 << " + " << nsteps << " * " << dt << " = " << (t_0 + nsteps * dt)
                               << " != " << t_end << ")");
        throw std::runtime_error("t_0 + nsteps * dt != t_end");
      }
    } else if (nsteps != 0) {
      t_end = t_0 + dt * nsteps;
    }
    const size_t niter = get_value<size_t>("num_iters", 10);

    //auto rueck = pfasst::examples::heat_FE::run_sdc(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter);


using pfasst::quadrature::quadrature_factory;



	double t1, t2; 
	t1 = MPI_Wtime();

        auto sdc = std::make_shared<heat_FE_sdc_t>();
        auto FinEl   = make_shared<fe_manager>(nelements,1); 
        auto sweeper = std::make_shared<sweeper_t>(FinEl, 0);


        sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);

	
        sdc->add_sweeper(sweeper);

        sdc->set_options();

        sdc->status()->time() = t_0;
        sdc->status()->dt() = dt;
        sdc->status()->t_end() = t_end;
        sdc->status()->max_iterations() = niter;

        sdc->setup();

        sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());

        sdc->run();

        sdc->post_run();


        std::cout << "####################################################################################################################  nach post run "  <<  std::endl;
        auto anfang    = sweeper->exact(0);
	        std::cout << "##########################################################################################################################nach exact "  <<  std::endl;
        auto naeherung = sweeper->get_end_state();
	        std::cout << "#########################################################################################################################nach neaherung "  <<  std::endl;
        auto exact     = sweeper->exact(t_end);
	        std::cout << "########################################################################################nach exact "  <<  std::endl;
        

        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
	if(my_rank==0) for (int i=0; i< nelements/num_pro; i++){
	  	        std::cout << "##################################################  in schleife "  <<  std::endl;
          std::cout << Dune::PDELab::Backend::native(anfang->data())[i] << " " << Dune::PDELab::Backend::native(naeherung->data())[i] << "   " << Dune::PDELab::Backend::native(exact->data())[i]<< " "  <<  std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
	/*if(my_rank==1) for (int i=0; i< nelements/2; i++){
	  	        std::cout << "##################################################  in schleife "  <<  std::endl;
          std::cout << Dune::PDELab::Backend::native(anfang->data())[i] << " " << Dune::PDELab::Backend::native(naeherung->data())[i] << "   " << Dune::PDELab::Backend::native(exact->data())[i]<< " "  <<  std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);*/

	std::cout << "##################################################  ausgabe zu ende "  <<  std::endl;
        //sweeper->states()[sweeper->get_states().size()-1]->scaled_add(-1.0 , sweeper->exact(t_end));
	
	
        //if (my_rank==0) std::cout << "Error " << sweeper->states()[sweeper->get_states().size()-1]->get_data().infinity_norm() <<  std::endl ;

	t2 = MPI_Wtime(); 
	printf( "Elapsed time is %f\n", t2 - t1 ); 
	return 0;

  }


