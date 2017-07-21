#include <config.h>
#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include "../../datatypes/dune_vec.hpp"
#include "fischer_sweeper.hpp"



using namespace pfasst::examples::heat_FE;

using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>; 


using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  
using sweeper_t = test<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>,   FE_function >;

using pfasst::transfer_traits;
using pfasst::contrib::SpectralTransfer;
using pfasst::SDC;
using pfasst::quadrature::QuadratureType;
using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

using pfasst::config::get_value;
using pfasst::quadrature::QuadratureType;
using pfasst::quadrature::quadrature_factory;

using pfasst::examples::heat_FE::test;



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


    auto sdc = std::make_shared<heat_FE_sdc_t>();
	
    auto FinEl   = make_shared<fe_manager>(nelements,1); 

    auto sweeper = std::make_shared<sweeper_t>(FinEl->get_basis(0), 0, FinEl->get_grid());
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

    sdc->run();
    //do not need a post run for GaussRadau nodes not sure about that should ask robert
    //sdc->post_run();

    //auto naeherung = sweeper->get_end_state()->data();
    //auto exact     = sweeper->exact(t_end)->data();

    // calculate the error of the simulation
    auto A = sweeper->get_A_dune(); A*=-1.0; //get the stiffnessmatrix, which is defined negative (need to change that!) 
    auto tmp = sweeper->get_end_state()->data();
    sweeper->get_end_state()->scaled_add(-1.0 , sweeper->exact(t_end));
    A.mv(sweeper->get_end_state()->data(), tmp);
    std::cout << "error in H1 norm: " << tmp*sweeper->get_end_state()->data() << std::endl;
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


