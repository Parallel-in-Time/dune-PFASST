#include <memory>
#include <iostream>

#include <vector>

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>

#include "FE_sweeper.hpp"

#include <pfasst/encap/dune_vec.hpp>


using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>; 

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
                                       const double& dt, const double& t_end, const size_t niter, double newton, bool output, double tol)
      {
        using pfasst::quadrature::quadrature_factory;

        auto sdc = std::make_shared<heat_FE_sdc_t>();

	auto FinEl   = make_shared<fe_manager>(nelements,1); //das Gitter soll ein mal vergroebert werden

        auto sweeper = std::make_shared<sweeper_t>(FinEl, 0); //sdc soll auf dem verfeinerten Gitter laufen (level 0)

        sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
	sweeper->newton=newton; //genauigkeit fuer den inneren Newton Loeser
        sweeper->set_abs_residual_tol(tol);

        sdc->add_sweeper(sweeper);

        sdc->set_options();

        sdc->status()->time() = t_0;
        sdc->status()->dt() = dt;
        sdc->status()->t_end() = t_end;
        sdc->status()->max_iterations() = niter;

        sdc->setup();

        sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());
	Dune::BlockVector<Dune::FieldVector<double, 1> > w = sweeper->initial_state()->data();
	
        sdc->run();

        sdc->post_run();



	if(output){
        	auto naeherung = sweeper->get_end_state()->data();
        	auto exact     = sweeper->exact(t_end)->data();

        	for (int i=0; i< sweeper->get_end_state()->data().size(); i++){
          		std::cout << sweeper->exact(0)->data()[i] << " " << naeherung[i] << "   " << exact[i] << std::endl;
        	}
        }

        sweeper->get_end_state()->scaled_add(-1.0 , sweeper->exact(t_end));
	std::cout << "ERROR (maximal difference of one component of the computed solution to analytic solution): " << sweeper->get_end_state()->norm0()<< std::endl;
	std::cout << "the corresponding linear system were solved " << sweeper->num_solves << " times" << std::endl; 
	std::cout << "(you solve this system in every time step for every time node for every outer iteration and for every Newton iteration)" << std::endl ;
	std::cout << "groesse loesungsvektor " << sweeper->get_end_state()->data().size() << std::endl ;
	std::cout << "Parameter " << sweeper->_n << " " << sweeper->_nu << std::endl ;

	
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

    pfasst::init(argc, argv, sweeper_t::init_opts);

    const size_t nelements = get_value<size_t>("num_elements", 180); //Anzahl der Elemente pro Dimension
    const size_t nnodes = get_value<size_t>("num_nodes", 3); 
    const QuadratureType quad_type = QuadratureType::GaussRadau;
    const double t_0 = 0.0;
    const double dt = get_value<double>("dt", 0.05);
    double t_end = get_value<double>("tend", 0.1);
    double newton = get_value<double>("newton", 0.1);
    bool output = get_value<double>("output", 0);
    const size_t niter = get_value<size_t>("num_iters", 10);
    double tol = get_value<double>("abs_res_tol", 1e-12);

    pfasst::examples::heat_FE::run_sdc(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter, newton, output, tol);

  }

#endif
