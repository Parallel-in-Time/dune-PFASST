//#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>


#include <config.h>

#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include <mpi.h>

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>

#include <pfasst/comm/mpi_p2p.hpp>
#include <pfasst/controller/two_level_pfasst.hpp>



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



#include "FE_sweeper.hpp"

#include "../../datatypes/dune_vec.hpp"
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
        typedef Grid::LeafGridView GV;
        //GV gv=gridp->leafGridView();
	//typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
	typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF,double, 1> FEM;  
        //FEM fem(gv);

  	// Make grid function space
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	//GFS gfs(gv,fem);

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double,  GFS, MatrixType>;
using pfasst::encap::DuneEncapsulation;


using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::contrib::SpectralTransfer;
using pfasst::TwoLevelPfasst;
typedef pfasst::comm::MpiP2P CommunicatorType;

using pfasst::examples::heat_FE::Heat_FE;

typedef DuneEncapsulation<double, double, GFS, MatrixType>                     EncapType;


typedef Heat_FE<pfasst::sweeper_traits<encap_traits_t>, MatrixType> SweeperType;


typedef pfasst::transfer_traits<SweeperType, SweeperType, 1>       TransferTraits;
typedef SpectralTransfer<TransferTraits>                           TransferType;


namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      void run_pfasst(const size_t nelements, const size_t basisorder, const size_t dim, const size_t& nnodes, const pfasst::quadrature::QuadratureType& quad_type,
                      const double& t_0, const double& dt, const double& t_end, const size_t& niter)
      {

        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

        TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
        pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);

	
        auto FinEl = make_shared<fe_manager>(nelements, 2);

        auto coarse = std::make_shared<SweeperType>(FinEl, 1);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<SweeperType>(FinEl, 0);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        coarse->is_coarse= true;
        fine->is_coarse=false;
        
        auto transfer = std::make_shared<TransferType>();


        transfer->create(FinEl);



	
	

        pfasst.add_transfer(transfer);
	
	
        pfasst.add_sweeper(coarse, true);
        pfasst.add_sweeper(fine);
        

        pfasst.set_options();


        pfasst.status()->time() = t_0;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_end;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

        coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());
        fine->initial_state() = fine->exact(pfasst.get_status()->get_time());

        pfasst.run();
        pfasst.post_run();



        MPI_Barrier(MPI_COMM_WORLD);


        
        if(my_rank==num_pro-1) {
        auto anfang    = fine->exact(0)->data();
        auto naeherung = fine->get_end_state()->data();
        auto exact     = fine->exact(t_end)->data();
        for (int i=0; i< fine->get_end_state()->data().size(); i++){
          std::cout << anfang[i] << " " << naeherung[i] << "   " << exact[i] << " "  <<  std::endl;
        }

        std::cout << "******************************************* " << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        fine->get_end_state()->scaled_add(-1.0, fine->exact(t_end));
        std::cout << "Fehler: "  << fine->get_end_state()->norm0() << " " << std::endl;


        std::cout << "******************************************* " << std::endl;

        
        
        ofstream f;
        stringstream ss;
        ss << nelements;
        string s = "solution_pfasst/" + ss.str() + ".dat";
        f.open(s, ios::app | std::ios::out );
        f << nelements << " " << dt << " "<< fine->get_end_state()->norm0() << endl;
        //f << nelements << " " << dt << " "<< x.infinity_norm()<< endl;

        f.close();

        ofstream ff;
        stringstream sss;
        sss << nelements << "_iter";
        string st = "solution_pfasst/" + sss.str() + ".dat";
        ff.open(st, ios::app | std::ios::out );
        auto iter = pfasst._it_per_step;
        for (const auto &line : iter) {
          ff << dt <<"  " << line << std::endl;
        }

        ff.close();
        
        
        
        }
      }
    }  // ::pfasst::examples::heat_FE
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  //Dune::MPIHelper::instance(argc, argv);
  MPI_Init(&argc, &argv);

  double t1, t2; 
  t1 = MPI_Wtime(); 

  pfasst::init(argc, argv, SweeperType::init_opts);
  pfasst::Status<double>::create_mpi_datatype();


  const size_t nelements = get_value<size_t>("num_elements", 10); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.1);
  double t_end = get_value<double>("tend", 0.8);
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
  const size_t niter = get_value<size_t>("num_iters", 50);

  pfasst::examples::heat_FE::run_pfasst(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter);

  pfasst::Status<double>::free_mpi_datatype();
  t2 = MPI_Wtime(); 
  printf( "Elapsed time is %f\n", t2 - t1 ); 
  MPI_Finalize();

  return 0;
}
