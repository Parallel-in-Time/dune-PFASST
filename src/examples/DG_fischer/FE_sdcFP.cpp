#include <config.h>
#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include "dune_vec_multi.hpp"
//#include "../../datatypes/dune_vec.hpp"
#include "fischer_sweeper.hpp"

//#include "assemble.hpp"

#include <dune/fufem/functionspacebases/p1nodalbasis.hh>
#include <dune/fufem/assemblers/operatorassembler.hh>
#include <dune/fufem/assemblers/functionalassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/l2functionalassembler.hh>
#include <dune/fufem/functiontools/basisinterpolator.hh>


using namespace pfasst::examples::fischer_example;

using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double,DIMENSION,blockSize>; 

//using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  
using FE_function = Dune::Functions::DGQkGLBlockBasis<GridType::LeafGridView, k>;

using sweeper_t = fischer_sweeper<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>,   FE_function >;

using pfasst::transfer_traits;
using pfasst::contrib::SpectralTransfer;
using pfasst::SDC;
using pfasst::quadrature::QuadratureType;
using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

using pfasst::config::get_value;
using pfasst::quadrature::QuadratureType;
using pfasst::quadrature::quadrature_factory;

using pfasst::examples::fischer_example::fischer_sweeper;





int main(int argc, char** argv) {
      
    Dune::MPIHelper::instance(argc, argv);
    
    
    // read in the parameter from terminal and make initial setup
    
    pfasst::init(argc, argv, sweeper_t::init_opts);

    const size_t nelements = get_value<size_t>("num_elements", 180);    // spacial dimension: number of grid points per dimension on the coase level
    
    const double t_0 = 0.0;                                             // left point of the time intervall is zero 
    const double dt = get_value<double>("dt", 0.05);                    // size of timesteping
    double t_end = get_value<double>("tend", 0.1);                      // right point of the time intervall  
    const size_t nnodes = get_value<size_t>("num_nodes", 3);            // time intervall: number of sdc quadrature points
    const QuadratureType quad_type = QuadratureType::GaussRadau;        // quadrature type
    const size_t niter = get_value<size_t>("num_iters", 10);            // maximal number of sdc iterations

    
    
     auto gridptr = Dune::StructuredGridFactory<GridType>::createCubeGrid({0,0},{1,1},{{2,2}}); //({0},{1},{{2}});
     gridptr->globalRefine(4);
     
     //std::shared_ptr<GridType> grid = std::make_shared<GridType>(gridptr);

     
     auto basis2 = FE_function{gridptr->leafGridView()};
//     
     std::shared_ptr<FE_function> basis = std::make_shared<FE_function>(basis2);
    
    //typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
    
    //typedef GridType::LevelGridView GridView;
    //using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
    
    std::shared_ptr<TransferOperatorAssembler<GridType>> transfer;

    std::shared_ptr<std::vector<MatrixType*>> transferMatrix;

//     std::shared_ptr<GridType> grid;
//     
// 
// //     int n_levels=2;
// // 
// //     std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); ; 
// // 
// //     
//      Dune::FieldVector<double,DIMENSION> hR = {200};
//      Dune::FieldVector<double,DIMENSION> hL = {-200};
//      array<int,DIMENSION> n;
//      std::fill(n.begin(), n.end(), nelements); 	    
//  #if HAVE_MPI
//      grid = std::make_shared<GridType>(hL, hR, n, std::bitset<DIMENSION>{0ULL}, 1, MPI_COMM_SELF);
//  #else
//      grid = std::make_shared<GridType>(hL, hR, n);
//  #endif
//     for (int i=0; i<n_levels; i++){	      
// 	      grid->globalRefine((bool) i);
// 	      auto view = grid->levelGridView(i);
//               fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); //grid->levelGridView(i));//gridView);
//     } 
    

    
    
    
    
    
    
    
    
    
    
//      typedef P1NodalBasis<GridView,double> FEBasis;
//      FEBasis feBasis(grid->levelGridView(1));
// 
//      // Build stiffness matrix
//      LaplaceAssembler<GridType, FEBasis::LocalFiniteElement, FEBasis::LocalFiniteElement> laplaceStiffness;
//      MatrixType stiffnessMatrix;
//      stiffnessMatrix*=-1;
//      
//      OperatorAssembler<FEBasis,FEBasis> operatorAssembler(feBasis, feBasis);
//      operatorAssembler.assemble(laplaceStiffness, stiffnessMatrix);
// 
//      MassAssembler<GridType, FEBasis::LocalFiniteElement, FEBasis::LocalFiniteElement> localMass;
//      MatrixType massMatrix;
//      operatorAssembler.assemble(localMass, massMatrix);
    
    
    //std::shared_ptr<FE_function> basis = std::make_shared<FE_function>(grid->leafGridView());

    
    
    
    
    
    auto sdc = std::make_shared<heat_FE_sdc_t>();
	

    
    MatrixType mass;
    MatrixType stiffness;

    //std::shared_ptr<BasisFunction> basis;   
    //basis =fe_basis2;
    //assembleProblem(basis, stiffness, mass);
    //std::cout << "vor " <<std::endl;
    
    auto sweeper = std::make_shared<sweeper_t>(basis , 0, gridptr); // mass and stiff are just dummies
    sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
    
    sdc->add_sweeper(sweeper);
    sdc->set_options();
    sdc->status()->time() = t_0;
    sdc->status()->dt() = dt;
    sdc->status()->t_end() = t_end;
    sdc->status()->max_iterations() = niter;
    sdc->setup();

    // initial value is given by the exact solution on time step 0
    sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());

    
    
    auto xBE = Dune::Fufem::istlVectorBackend(sweeper->initial_state()->data());
    auto xFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>( FE_function{gridptr->leafGridView()}, Dune::TypeTree::hybridTreePath(), xBE);
    Dune::VTKWriter<typename GridType::LeafGridView> vtkWriter(gridptr->leafGridView());

    vtkWriter.addVertexData(xFunction, Dune::VTK::FieldInfo("x", Dune::VTK::FieldInfo::Type::scalar, 1));
    const std::string filename="initial";
    vtkWriter.write(filename);
    
    for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << sweeper->initial_state()->data()[i] << std::endl;
    //std::exit(0);
    sdc->run();
    sdc->post_run();

    
    auto xBE2 = Dune::Fufem::istlVectorBackend(sweeper->get_end_state()->data());
    auto xFunction2 = Dune::Functions::makeDiscreteGlobalBasisFunction<double>( FE_function{gridptr->leafGridView()}, Dune::TypeTree::hybridTreePath(), xBE2);
    Dune::VTKWriter<typename GridType::LeafGridView> vtkWriter2(gridptr->leafGridView());

    vtkWriter2.addVertexData(xFunction2, Dune::VTK::FieldInfo("x", Dune::VTK::FieldInfo::Type::scalar, 1));
    const std::string filename2="end";
    vtkWriter2.write(filename2);
    
    for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << sweeper->get_end_state()->data()[i] <<  std::endl;

    std::cout << "error in infinity norm: " << sweeper->get_end_state()->get_data().infinity_norm() <<  std::endl ;

   

 
}


