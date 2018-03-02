
#include <memory>
#include <iostream>

#include <vector>

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>


#include "FE_sweeper.hpp"
#include "dune_vec_multi.hpp"
#include "../../finite_element_stuff/spectral_transfer_gs.hpp"

//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

const size_t DIM = 2;            //Räumliche Dimension des Rechengebiets ruth_dim

const size_t BASIS_ORDER = 1;    //maximale Ordnung der Lagrange Basisfunktionen

const size_t NR_OF_COMP = 2;
//////////////////////////////////////////////////////////////////////////////////////

using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, DIM, NR_OF_COMP>; //ruth dim, NR_OF_COMP



 namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      using sweeper_t = Heat_FE<dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM, NR_OF_COMP>>;
      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::SDC;
      using pfasst::quadrature::QuadratureType;
      using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;  //?2

      shared_ptr<heat_FE_sdc_t> run_sdc(const size_t nelements, const size_t basisorder, const size_t dim, const size_t nnodes,
                                       const QuadratureType& quad_type, const double& t_0,
                                       const double& dt, const double& t_end, const size_t niter)
      {
        using pfasst::quadrature::quadrature_factory;

        auto sdc = std::make_shared<heat_FE_sdc_t>();

        //auto sweeper = std::make_shared<sweeper_t>(nelements, basisorder, 0);



        auto FinEl   = make_shared<fe_manager>(nelements,2);
        auto sweeper = std::make_shared<sweeper_t>(FinEl, 1);

        sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        sdc->add_sweeper(sweeper);

        sdc->set_options();

        sdc->status()->time() = t_0;
        sdc->status()->dt() = dt;
        sdc->status()->t_end() = t_end;
        sdc->status()->max_iterations() = niter;

        sdc->setup();

        sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());

	Dune::BlockVector<Dune::FieldVector<double, NR_OF_COMP> > w = sweeper->initial_state()->data();
	
 	double time1=0.0, tstart;      // time measurment variables
 
 	tstart = clock();              // start

        sdc->run();

        sdc->post_run();

 	time1 += clock() - tstart;     // end..
 
 	time1 = time1/CLOCKS_PER_SEC;  // rescale to seconds

 	cout << "  time = " << time1 << " sec." << endl;
	double t1, t2;
	t1=MPI_Wtime();	



        ofstream ff;
        stringstream sss;
        sss << nelements << "_" << dt << "_" << t_end;
        string st = "solution_sdc/" + sss.str() + ".dat";
        ff.open(st, ios::app | std::ios::out );
        auto iter = sdc->_it_per_step;
	int i=0;
        for (const auto &line : iter) {
          ff << i <<"  " << line << std::endl;
          i++;
        }
 	ff << "  time = " << time1 << endl;

        ff << "-----------------------------------------------  " << std::endl;
        ff.close();

	
        //return sdc;
      }
    }
  }
}






#ifndef PFASST_UNIT_TESTING
  int main(int argc, char** argv) {
    using pfasst::config::get_value;
    using pfasst::quadrature::QuadratureType;
    using pfasst::examples::heat_FE::Heat_FE;

    using sweeper_t = Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM, NR_OF_COMP>>;

    pfasst::init(argc, argv, sweeper_t::init_opts);

    const size_t nelements = get_value<size_t>("num_elements", 200); //Anzahl der Elemente pro Dimension
    const size_t nnodes = get_value<size_t>("num_nodes",3);
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
    
    std::cout << "nsteps = " << nsteps << std::endl; 
    
    const size_t niter = get_value<size_t>("num_iters", 100);

    pfasst::examples::heat_FE::run_sdc(nelements, BASIS_ORDER, DIM, nnodes, quad_type, t_0, dt, t_end, niter);
    

 
  return 0;
    
    
 
  }

#endif
