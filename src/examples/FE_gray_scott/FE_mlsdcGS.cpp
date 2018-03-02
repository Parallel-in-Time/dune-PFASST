
#include <fenv.h>
#include <memory>
#include <stdexcept>
#include <vector>

using std::shared_ptr;

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/two_level_mlsdc.hpp>

#include "FE_sweeper.hpp"
#include "dune_vec_multi.hpp"
#include "../../finite_element_stuff/spectral_transfer_gs.hpp"










//using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 2, 2>;


// using encap_traits_t = pfasst::eigen3_encap_traits<double, double, 2>;


//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

const size_t DIM = 2;            //Räumliche Dimension des Rechengebiets ruth_dim

const size_t BASIS_ORDER = 1;    //maximale Ordnung der Lagrange Basisfunktionen

const size_t NR_OF_COMP=2;
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

      //using sweeper_t = Heat_FE<pfasst::sweeper_traits<encap_traits_t, BASIS_ORDER, DIM, NR_OF_COMP>>;
      using sweeper_t = Heat_FE<   dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM, NR_OF_COMP>>;
      using transfer_traits_t = pfasst::transfer_traits<sweeper_t, sweeper_t, 1>;
      using transfer_t = SpectralTransfer<transfer_traits_t>;
      using heat_FE_mlsdc_t = TwoLevelMLSDC<transfer_t>;

      //shared_ptr<heat_FE_mlsdc_t>
      void run_mlsdc(const size_t nelements, const size_t basisorder, const size_t DIM, const size_t coarse_factor,
                                           const size_t nnodes, const QuadratureType& quad_type,
                                           const double& t_0, const double& dt, const double& t_end,
                                           const size_t niter) {

        auto mlsdc = std::make_shared<heat_FE_mlsdc_t>();

        auto FinEl = make_shared<fe_manager>(nelements, 2); //spatial_dofs

        using pfasst::quadrature::quadrature_factory;

        auto coarse = std::make_shared<sweeper_t>(FinEl, 1);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<sweeper_t>(FinEl, 0);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);
	
        std::cout <<"transfer" << std::endl;
        auto transfer = std::make_shared<transfer_t>();
        transfer->create(FinEl);
        //mlsdc->add_sweeper(coarse, true);
        //mlsdc->add_sweeper(fine, false);


	std::cout <<"vor add_transfer" << std::endl;
        mlsdc->add_sweeper(coarse, true);
        mlsdc->add_sweeper(fine);
        mlsdc->add_transfer(transfer);

        std::cout <<"vor set_options" << std::endl;
        mlsdc->set_options();


        mlsdc->status()->time() = t_0;
        mlsdc->status()->dt() = dt;
        mlsdc->status()->t_end() = t_end;
        mlsdc->status()->max_iterations() = niter;

	std::cout <<"vor setup" << std::endl;
        mlsdc->setup();

        coarse->is_coarse= true;
        fine->is_coarse=false;

        coarse->initial_state() = coarse->exact(mlsdc->get_status()->get_time());
        fine->initial_state() = fine->exact(mlsdc->get_status()->get_time());

        /*for (int i=0; i< fine->initial_state()->data().size(); i++){
          std::cout << "Anfangswerte feiner Sweeper: " << " " << fine->initial_state()->data()[i] << std::endl;
        }

        std::cout  <<  std::endl;

        for (int i=0; i< coarse->initial_state()->data().size(); i++){
          std::cout << "Anfangswerte grober Sweeper: " << " " << coarse->initial_state()->data()[i] <<  std::endl;
        }*/


 	double time1=0.0, tstart;      // time measurment variables
 
 	tstart = clock();              // start

	     
        //std::cout << "*********************************vor run"<<  std::endl ;
        mlsdc->run();
        //std::cout << "*********************************nach run"<<  std::endl ;

        mlsdc->post_run();


 	time1 += clock() - tstart;     // end

 	time1 = time1/CLOCKS_PER_SEC;  // rescale to seconds

 	cout << "  time = " << time1 << " sec." << endl;

        ofstream ff;
        stringstream sss;
        sss << nelements*2 << "_" << dt << "_" << t_end;
        string st = "solution_mlsdc/" + sss.str() + ".dat";
        ff.open(st, ios::app | std::ios::out );
        auto iter = mlsdc->_it_per_step;
	int i=0;
        for (const auto &line : iter) {
          ff << i <<"  " << line << std::endl;
          i++;
        }
 	ff << "  time = " << time1 << endl;

        ff << "-----------------------------------------------  " << std::endl;
        ff.close();
        

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
  //using sweeper_t = pfasst::examples::heat_FE::Heat_FE<pfasst::sweeper_traits<encap_traits_t, BASIS_ORDER, DIM, NR_OF_COMP>>;
  using sweeper_t = pfasst::examples::heat_FE::Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM, NR_OF_COMP>>;

  pfasst::init(argc, argv, sweeper_t::init_opts);

  const size_t nelements = get_value<size_t>("num_elements", 32); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  //const size_t ndofs = get_value<size_t>("num_dofs", 8);
  const size_t coarse_factor = get_value<size_t>("coarse_factor", 1);
  //const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 2);
  double t_end = get_value<double>("tend", 6);
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
  const size_t niter = get_value<size_t>("num_iters", 20);

  pfasst::examples::heat_FE::run_mlsdc(nelements, BASIS_ORDER, DIM, coarse_factor, nnodes, quad_type, t_0, dt, t_end, niter);
}
#endif 
