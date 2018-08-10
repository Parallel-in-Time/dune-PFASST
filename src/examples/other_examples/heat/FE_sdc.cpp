#include "config.h"

#include <memory>
#include <iostream>
#include <vector>

//#include "dune_includes"
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
//create mass + dt stiffness 
#include"nonlinearpoissonfem.hh"
//creates mass
#include"mass_operator.hh"
//constants


#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/config.hpp>

#include "FE_sweeper.hpp"
#include "spectral_transfer.hpp"





//Dune::PDELab::Backend::Vector<GFS,SpatialPrecision>

typedef double RF; 
typedef Dune::YaspGrid<1> GridType;
typedef GridType::ctype DF;
typedef GridType::LeafGridView GV;
typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
typedef Dune::PDELab::istl::VectorBackend<> VBE;
typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

typedef typename GFS::template
ConstraintsContainer<RF>::Type CC;
typedef NonlinearPoissonFEM<FEM> LOP;
typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GO;
typedef massFEM<FEM> LOPm;
typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOPm,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GOm;
typedef GOm MatrixType;

using std::shared_ptr;
using dune_vector_encap_traits_t = pfasst::encap::dune_vec_encap_traits<GFS, double, double, 1>;
using namespace pfasst::examples::heat_FE;
using pfasst::config::get_value;
using pfasst::quadrature::QuadratureType;
using pfasst::examples::heat_FE::Heat_FE;
using sweeper_t = Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<dune_vector_encap_traits_t, 1, 1>, MatrixType>;
using pfasst::transfer_traits;
using pfasst::contrib::SpectralTransfer;
using pfasst::SDC;
using pfasst::quadrature::QuadratureType;
using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;
using pfasst::quadrature::quadrature_factory;

int main(int argc, char** argv) {


    pfasst::init(argc, argv, sweeper_t::init_opts);

    auto const   nelements = pfasst::config::get_value<size_t>("num_elements", 3); 
    const size_t nnodes    = get_value<size_t>("num_nodes", 3);
    const QuadratureType quad_type = QuadratureType::GaussRadau;
    const double t_0 = 0.0;
    const double dt = get_value<double>("dt", 0.05);
    double t_end = get_value<double>("tend", 0.1);
    const size_t niter = get_value<size_t>("num_iters", 10);


    auto sdc = std::make_shared<heat_FE_sdc_t>();
    //auto FinEl   = make_shared<fe_manager>(nelements,1); 
    auto sweeper = std::make_shared<sweeper_t>(0); //make a coarse sweeper

    sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
	
    sdc->add_sweeper(sweeper);

    sdc->set_options();

    sdc->status()->time() = t_0;
    sdc->status()->dt() = dt;
    sdc->status()->t_end() = t_end;
    sdc->status()->max_iterations() = niter;

    sdc->setup();

    //sweeper->initial_state() = 
    sweeper->exact(sdc->get_status()->get_time());

    sdc->run();

    sdc->post_run();



    //sweeper->states()[sweeper->get_states().size()-1]->scaled_add(-1.0 , sweeper->exact(t_end));
	
	
    //std::cout << "Error " << sweeper->states()[sweeper->get_states().size()-1]->get_data().infinity_norm() <<  std::endl ;

}








