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
#include <pfasst/controller/two_level_mlsdc.hpp>

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
        //typedef Grid::LeafGridView GV;
	typedef Grid::LevelGridView GV;
        //GV gv=gridp->leafGridView();
	typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM;
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
      using pfasst::TwoLevelMLSDC;
      using pfasst::quadrature::QuadratureType;
      using transfer_traits_t = pfasst::transfer_traits<sweeper_t, sweeper_t, 1>;
      using transfer_t = SpectralTransfer<transfer_traits_t>;


      using heat_FE_mlsdc_t = TwoLevelMLSDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

 
    }
  }
}


int main(int argc, char** argv)
{

  //feenableexcept(FE_INVALID | FE_OVERFLOW);



  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;
  using sweeper_t = pfasst::examples::heat_FE::Heat_FE<pfasst::sweeper_traits<encap_traits_t>, MatrixType>;

	Dune::MPIHelper&
      	helper = Dune::MPIHelper::instance(argc, argv);

  //Dune::MPIHelper::instance(argc, argv);
  pfasst::init(argc, argv, sweeper_t::init_opts);

  const size_t nelements = get_value<size_t>("num_elements", 20);
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const size_t coarse_factor = get_value<size_t>("coarse_factor", 1);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.001);
  double t_end = get_value<double>("tend", 0.001);
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
        auto mlsdc = std::make_shared<heat_FE_mlsdc_t>();

        auto FinEl = make_shared<fe_manager>(nelements, 2); 



        using pfasst::quadrature::quadrature_factory;

        auto coarse = std::make_shared<sweeper_t>(FinEl, 0); //0
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<sweeper_t>(FinEl, 1);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        coarse->is_coarse= true;
        fine->is_coarse=false;

        auto transfer = std::make_shared<transfer_t>();


	const int dim=1;
	const int degree =1;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<double,dim> h = {1};	      
	std::array<int,dim> n;
	std::fill(n.begin(), n.end(), nelements);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(h,n));

        gridp->refineOptions(false); 
        gridp->globalRefine(1);
	typedef Grid::LevelGridView GV;
        GV gv1=gridp->levelGridView(1);
        GV gv0=gridp->levelGridView(0);
	//typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
	typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF,double, 1> FEM;  
        FEM fem0(gv0);
        FEM fem1(gv1);
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	GFS gfs0(gv0,fem0);
  	GFS gfs1(gv1,fem0);




        transfer->create(FinEl, gfs0, gfs1);




	
        mlsdc->add_sweeper(coarse, true);
        mlsdc->add_sweeper(fine);

        mlsdc->add_transfer(transfer);

	
        mlsdc->set_options();


        mlsdc->status()->time() = t_0;
        mlsdc->status()->dt() = dt;
        mlsdc->status()->t_end() = t_end;
        mlsdc->status()->max_iterations() = niter;


        mlsdc->setup();


        coarse->initial_state() = coarse->exact(mlsdc->get_status()->get_time());
        fine->initial_state() = fine->exact(mlsdc->get_status()->get_time());
       int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

	if(my_rank==0) for (int i=0; i< Dune::PDELab::Backend::native(coarse->initial_state()->data()).N(); i++){
          std::cout << "initial state " << Dune::PDELab::Backend::native(coarse->initial_state()->data())[i] <<  std::endl;
        }        MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==1) for (int i=0; i< Dune::PDELab::Backend::native(coarse->initial_state()->data()).N(); i++){
          std::cout << "initial state " << Dune::PDELab::Backend::native(coarse->initial_state()->data())[i] <<  std::endl;
        }        MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==0) for (int i=0; i< Dune::PDELab::Backend::native(fine->initial_state()->data()).N(); i++){
          std::cout << "initial state " << Dune::PDELab::Backend::native(fine->initial_state()->data())[i] <<  std::endl;
        }        MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==1) for (int i=0; i< Dune::PDELab::Backend::native(fine->initial_state()->data()).N(); i++){
          std::cout << "initial state " << Dune::PDELab::Backend::native(fine->initial_state()->data())[i] <<  std::endl;
        }        MPI_Barrier(MPI_COMM_WORLD);
	//std::exit(0);

        /*for (int i=0; i< fine->initial_state()->data().size(); i++){
          std::cout << "Anfangswerte feiner Sweeper: " << " " << fine->initial_state()->data()[i] << std::endl;
        }

        std::cout  <<  std::endl;

        for (int i=0; i< coarse->initial_state()->data().size(); i++){
          std::cout << "Anfangswerte grober Sweeper: " << " " << coarse->initial_state()->data()[i] <<  std::endl;
        }*/




        mlsdc->run();


        mlsdc->post_run();


 	std::cout << "####################################################################################################################  nach post run "  <<  std::endl;
        auto anfang    = fine->exact(0);
	        std::cout << "##########################################################################################################################nach exact "  <<  std::endl;
        auto naeherung = fine->get_end_state();
	        std::cout << "#########################################################################################################################nach neaherung "  <<  std::endl;
        auto exact     = fine->exact(t_end);
	        std::cout << "########################################################################################nach exact "  <<  std::endl;
        

 

	if(my_rank==0) for (int i=0; i< Dune::PDELab::Backend::native(anfang->data()).N(); i++){
          std::cout << "ergebnis 1 " << Dune::PDELab::Backend::native(anfang->data())[i] << " " << Dune::PDELab::Backend::native(naeherung->data())[i] << "   " << Dune::PDELab::Backend::native(exact->data())[i]<< " "  <<  std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==1) for (int i=0; i< nelements; i++){
          std::cout << "ergebnis 1 " << Dune::PDELab::Backend::native(anfang->data())[i] << " " << Dune::PDELab::Backend::native(naeherung->data())[i] << "   " << Dune::PDELab::Backend::native(exact->data())[i]<< " "  <<  std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        std::cout << "******************************************* " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << "Fehler: " <<  std::endl ;
        fine->states()[fine->get_states().size()-1]->scaled_add(-1.0 , fine->exact(t_end));
	double max = fine->states()[fine->get_states().size()-1]->norm0();
	double norm;
	MPI_Reduce(
	    &max,
	    &norm,
	    1,
	    MPI_DOUBLE,
	    MPI_MAX,
	    0,
	    MPI_COMM_WORLD);
        if(my_rank==0) std::cout <<"FEHLER" <<  max <<  std::endl ;
        std::cout << "******************************************* " <<  std::endl ;



        std::cout << "******************************************* " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << "Fehler: " <<  std::endl ;
        coarse->states()[coarse->get_states().size()-1]->scaled_add(-1.0 , coarse->exact(t_end));
	double cmax = coarse->states()[coarse->get_states().size()-1]->norm0();
        std::cout <<my_rank  <<"** norm ***************************************** " << cmax << std::endl ;
	double cnorm;
	MPI_Reduce(
	    &cmax,
	    &cnorm,
	    1,
	    MPI_DOUBLE,
	    MPI_MAX,
	    0,
	    MPI_COMM_WORLD);
        if(my_rank==0) std::cout << cnorm <<  std::endl ;
        std::cout << "******************************************* " <<  std::endl ;



}



