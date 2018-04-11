#include <config.h>
#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include <dune/istl/schwarz.hh> // Dune::OverlappingSchwarzScalarProduct, Dune::OverlappingSchwarzOperator


#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication


#include "../../datatypes/dune_vec.hpp"
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

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>; 


using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  
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

    
    
    typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
    typedef GridType::LevelGridView GridView;
    using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
    
    std::shared_ptr<TransferOperatorAssembler<GridType>> transfer;

    std::shared_ptr<std::vector<MatrixType*>> transferMatrix;

    std::shared_ptr<GridType> grid;
    std::cout << "test "<< std::endl;

    int n_levels=2;

    std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); ; 

    
    Dune::FieldVector<double,DIMENSION> hR = {200};
    Dune::FieldVector<double,DIMENSION> hL = {-200};
    array<int,DIMENSION> n;
    std::fill(n.begin(), n.end(), nelements); 	    
#if HAVE_MPI
    grid = std::make_shared<GridType>(hL, hR, n, std::bitset<DIMENSION>{0ULL}, 1, MPI_COMM_SELF);
#else
    grid = std::make_shared<GridType>(hL, hR, n);
#endif
    for (int i=0; i<n_levels; i++){	      
	      grid->globalRefine((bool) i);
	      auto view = grid->levelGridView(i);
              fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); //grid->levelGridView(i));//gridView);
    } 
    

    
    
    
    
    
    
    
    
    
    
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
    
    
    
    
    
    
    
    
    auto sdc = std::make_shared<heat_FE_sdc_t>();
	

    
    MatrixType mass;
    MatrixType stiffness;

    //std::shared_ptr<BasisFunction> basis;   
    //basis =fe_basis2;
    //assembleProblem(basis, stiffness, mass);
    //std::cout << "vor " <<std::endl;
    
    auto sweeper = std::make_shared<sweeper_t>(fe_basis[0] , 0, grid); // mass and stiff are just dummies
    sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
    
    sdc->add_sweeper(sweeper);
    sdc->set_options();
    sdc->status()->time() = t_0;
    sdc->status()->dt() = dt;
    sdc->status()->t_end() = t_end;
    sdc->status()->max_iterations() = niter;
    sdc->setup();


		using DuneCommunication = Dune::OwnerOverlapCopyCommunication<std::size_t>; 
    // initial value is given by the exact solution on time step 0
    sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());

    for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << sweeper->initial_state()->data()[i] << std::endl;
    //std::exit(0);
    std::cout << "vor run "<< std::endl;
    sdc->run();
    std::cout << "nach run"<< std::endl;    
    
    //do not need a post run for GaussRadau nodes not sure about that should ask robert
    //sdc->post_run();

    auto naeherung = sweeper->get_end_state()->data();
    auto exact     = sweeper->exact(t_end)->data();
    auto initial   = sweeper->exact(0)->data();
    for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << initial[i] << " " << naeherung[i] << " " << exact[i] << std::endl;

    // calculate the error of the simulation
    //auto A = sweeper->get_A_dune(); A*=-1.0; //get the stiffnessmatrix, which is defined negative (need to change that!) 
    //auto tmp = sweeper->get_end_state()->data();
    sweeper->get_end_state()->scaled_add(-1.0 , sweeper->exact(t_end));
    //A.mv(sweeper->get_end_state()->data(), tmp);
    //std::cout << "error in H1 norm: " << tmp*sweeper->get_end_state()->data() << std::endl;
    std::cout << "error in infinity norm: " << sweeper->get_end_state()->norm0()<<  std::endl ;

    bool output=false;
    if (output){
        //write in a file the error
        ofstream f;
        stringstream ss;
        ss << nelements;
        string s = "solution_sdc/" + ss.str() + ".dat";
        f.open(s, ios::app | std::ios::out );
        f << nelements << " " << dt << " "<< sweeper->get_end_state()->norm0()<< endl;
        f.close();

        //write in a file the number of sdc iterations
        ofstream ff;
        stringstream sss;
        sss << nelements <<  "_iter";
        string st = "solution_sdc/" + sss.str() + ".dat";
        ff.open(st, ios::app | std::ios::out );
        auto iter = sdc->_it_per_step;
        for (const auto &line : iter) {
            ff << dt <<"  " << line << std::endl;
        }
        ff.close();
    }	


 
}


