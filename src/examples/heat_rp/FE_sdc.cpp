#include <config.h>
#include <memory>
#include <iostream>
#include <vector>
#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>

#include <pfasst/config.hpp>


//#include <dune/istl/schwarz.hh> // Dune::OverlappingSchwarzScalarProduct, Dune::OverlappingSchwarzOperator


//#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
//#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication



#include "../../datatypes/dune_vec.hpp"

#include "FE_sweeper.hpp"

#include "spectral_transfer.hpp"
using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;


const size_t nelements = 10;

namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      using sweeper_t = Heat_FE<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>>;
      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::SDC;
      using pfasst::quadrature::QuadratureType;
      using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

      shared_ptr<heat_FE_sdc_t> run_sdc(const size_t nelements, const size_t basisorder, const size_t dim, const size_t nnodes,
                                       const QuadratureType& quad_type, const double& t_0,
                                       const double& dt, const double& t_end, const size_t niter)
      {
        using pfasst::quadrature::quadrature_factory;

        auto sdc = std::make_shared<heat_FE_sdc_t>();
        auto FinEl   = make_shared<fe_manager>(nelements,1); 
        auto sweeper = std::make_shared<sweeper_t>(FinEl, 0);


        sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        //sweeper->set_abs_residual_tol(1e-8);
	
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




        sweeper->states()[sweeper->get_states().size()-1]->scaled_add(-1.0 , sweeper->exact(t_end));
	
	
        std::cout << "Error " << sweeper->states()[sweeper->get_states().size()-1]->get_data().infinity_norm() <<  std::endl ;


	  ofstream f;
	  stringstream ss;
	  ss << nelements;
	  string s = "solution_sdc/" + ss.str() + ".dat";
	  f.open(s, ios::app | std::ios::out );
	  f << nelements << " " << dt << " "<< sweeper->states()[sweeper->get_states().size()-1]->norm0()<< endl;


	  f.close();

        ofstream ff;
        stringstream sss;
        sss << nelements << "_iter";
        string st = "solution_sdc/" + sss.str() + ".dat";
        ff.open(st, ios::app | std::ios::out );
        auto iter = sdc->_it_per_step;
        for (const auto &line : iter) {
          ff << dt <<"  " << line << std::endl;
        }

        ff.close();



        return sdc;
      }
    }
  }
}



#ifndef PFASST_UNIT_TESTING
  int main(int argc, char** argv) {
    using pfasst::config::get_value;
    using pfasst::quadrature::QuadratureType;
    using pfasst::examples::heat_FE::Heat_FE;

    using sweeper_t = Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>>;

  	MPI_Init(&argc, &argv);

        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

		std::cout << "soviele anfang"<< my_rank << " von " << num_pro << std::endl; //<< u_seq[i] << std::endl;

         Dune::FakeMPIHelper::instance(argc, argv);


    pfasst::init(argc, argv, sweeper_t::init_opts);

    auto const   nelements = pfasst::config::get_value<size_t>("num_elements", 10); 
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

//		using DuneCommunication = Dune::OwnerOverlapCopyCommunication<std::size_t>;
		std::cout << "soviele "<< my_rank << " von " << num_pro << std::endl; //<< u_seq[i] << std::endl;


double starttime, endtime;
       starttime = MPI_Wtime();
    pfasst::examples::heat_FE::run_sdc(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter);
endtime   = MPI_Wtime();
       printf("That took %f seconds\n",endtime-starttime);

MPI_Finalize();
  }

#endif
