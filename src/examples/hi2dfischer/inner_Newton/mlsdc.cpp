
#include <fenv.h>

#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/two_level_mlsdc.hpp>




#include <vector>


#include "FE_sweeper.hpp"
#include <pfasst/encap/dune_vec.hpp>
#include "../2d_transfer/spectral_transfer.hpp"


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;



//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

const size_t DIM = 1;            //R??umliche Dimension des Rechengebiets ruth_dim

const size_t BASIS_ORDER = 1;    //maximale Ordnung der Lagrange Basisfunktionen

//////////////////////////////////////////////////////////////////////////////////////
const size_t nelements = 100;




namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::TwoLevelMLSDC;
      using pfasst::quadrature::QuadratureType;

      using sweeper_t_coarse = Heat_FE<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>>;
      using sweeper_t_fine = Heat_FE<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>>;
      //using sweeper_t = Heat_FE<pfasst::sweeper_traits<encap_traits_t>>;
      using transfer_traits_t = pfasst::transfer_traits<sweeper_t_coarse, sweeper_t_fine, 1>;
      using transfer_t = SpectralTransfer<transfer_traits_t>;
      using heat_FE_mlsdc_t = TwoLevelMLSDC<transfer_t>;


      void run_mlsdc(const size_t nelements, const size_t basisorder, const size_t coarse_factor,
                                           const size_t nnodes, const QuadratureType& quad_type,
                                           const double& t_0, const double& dt, const double& t_end,
                                           const size_t niter, const double newton, bool output) {
        auto mlsdc = std::make_shared<heat_FE_mlsdc_t>();

        auto FinEl = make_shared<fe_manager>(nelements, 2);
        //mlsdc->grid_builder(nelements);

        using pfasst::quadrature::quadrature_factory;

        auto coarse = std::make_shared<sweeper_t_coarse>(FinEl, 1);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<sweeper_t_coarse>(FinEl, 0);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        
        coarse->is_coarse=true;
        coarse->newton=newton;
        fine->is_coarse=false;
        fine->newton=newton;
        fine->output=output;
        fine->output_level=1;
        //coarse->is_coarse=true;
        //fine->is_coarse=false;
        
        
        auto transfer = std::make_shared<transfer_t>();
        transfer->create(FinEl);
	
        //mlsdc->add_sweeper(coarse, true);
        //mlsdc->add_sweeper(fine, false);

    
           //fine->set_abs_residual_tol(1e-12);
           //coarse->set_abs_residual_tol(1e-12);
    
    
           
        
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

        mlsdc->run();

        mlsdc->post_run();


	std::cout << "the corresponding linear system were solved " << fine->num_solves << " times" << std::endl; 
	std::cout << "(you solve this system in every time step for every time node for every outer iteration and for every Newton iteration)" << std::endl ;
	std::cout << "Groesse des Loesungsvektors: " << fine->get_end_state()->data().size() << std::endl ;

        //return mlsdc;

        
        /*std::cout <<  "fein" << std::endl;
        auto naeherung = fine->get_end_state()->data();
        auto exact     = fine->exact(t_end)->data();
        for (int i=0; i< fine->get_end_state()->data().size(); i++){
          std::cout << fine->exact(0)->data()[i] << " " << naeherung[i] << "   " << exact[i] << std::endl;
        }

        std::cout << "******************************************* " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << "Fehler: " <<  std::endl ;
        //auto norm =  fine->exact(t_end))->data();
        fine->states()[fine->get_states().size()-1]->scaled_add(-1.0 , fine->exact(t_end));
        std::cout << fine->states()[fine->get_states().size()-1]->norm0()<<  std::endl ;
        std::cout << "number states " << fine->get_states().size() << std::endl ;
        std::cout << "******************************************* " <<  std::endl ;

        //std::cout << "Fehler: " <<  std::endl ;
        //coarse->states()[coarse->get_states().size()-1]->scaled_add(-1.0 , coarse->exact(t_end));
        //std::cout << coarse->states()[coarse->get_states().size()-1]->norm0()<<  std::endl ;
	*/





      }

    }  // ::pfasst::examples::heat_FE
  } // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{

  feenableexcept(FE_INVALID | FE_OVERFLOW);

  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;
  //using sweeper_t      = pfasst::examples::heat_FE::Heat_FE<pfasst::sweeper_traits<encap_traits_t,1,1>>;
  using sweeper_t_fine = pfasst::examples::heat_FE::Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, 2, DIMENSION>>;

  pfasst::init(argc, argv, sweeper_t_fine::init_opts);

  const size_t nelements = get_value<size_t>("num_elements", 180); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  //const size_t ndofs = get_value<size_t>("num_dofs", 8);
  const size_t coarse_factor = get_value<size_t>("coarse_factor", 1);
  //const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.05);
  double t_end = get_value<double>("tend", 0.1);
  size_t nsteps = get_value<size_t>("num_steps", 0);
  double newton = get_value<size_t>("newton", 1e-6);
  bool output = get_value<double>("output", 0);
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
  
    MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
    	auto st = MPI_Wtime();
  pfasst::examples::heat_FE::run_mlsdc(nelements, BASIS_ORDER, coarse_factor, nnodes, quad_type, t_0, dt, t_end, niter, newton, output);
  
      	auto ut = MPI_Wtime()-st;
        double time;
        MPI_Allreduce(&ut, &time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        std::cout << "Zeit ist am ende " << time << std::endl;
    MPI_Finalize();
}
#endif 
