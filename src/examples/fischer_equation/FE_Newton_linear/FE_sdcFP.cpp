#include <memory>
#include <iostream>

#include <vector>

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
//#include <pfasst/encap/dune_vec.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>

#include "FE_sweeper.hpp"

#include "../../datatypes/dune_vec.hpp"

//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

const size_t DIM = 1;            //Raeumliche Dimension des Rechengebiets ruth_dim

const size_t BASIS_ORDER = 1;    //maximale Ordnung der Lagrange Basisfunktionen

//////////////////////////////////////////////////////////////////////////////////////


using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>; 

namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      using sweeper_t = Heat_FE<dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM>>;
      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::SDC;
      using pfasst::quadrature::QuadratureType;
      using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

      shared_ptr<heat_FE_sdc_t> run_sdc(const size_t nelements, const size_t basisorder, const size_t dim, const size_t nnodes,
                                       const QuadratureType& quad_type, const double& t_0,
                                       const double& dt, const double& t_end, const size_t niter, double newton)
      {
        using pfasst::quadrature::quadrature_factory;

        auto sdc = std::make_shared<heat_FE_sdc_t>();
	
	auto FinEl   = make_shared<fe_manager>(nelements,1); 

        auto sweeper = std::make_shared<sweeper_t>(FinEl, 0);

        sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
	sweeper->newton=newton;

        //sweeper->set_abs_residual_tol(1e-10);
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


        /*if(BASIS_ORDER==1) {
          auto grid = (*sweeper).get_grid();
          typedef GridType::LeafGridView GridView;
          GridType::LeafGridView gridView = grid->leafGridView();
          VTKWriter<GridView> vtkWriter(gridView);
          typedef Dune::BlockVector<Dune::FieldVector<double, 1> > VectorType;
          VectorType x = sweeper->get_end_state()->data();
	  VectorType y = sweeper->exact(t_end)->data();
	  VectorType z = sweeper->initial_state()->data();
	 
	  
	    
          vtkWriter.addVertexData(x, "fe_solution");
          vtkWriter.addVertexData(y, "exact_solution");
          vtkWriter.addVertexData(w, "initial_data");

          vtkWriter.write("flame_result");
        }*/
	
	
	
	/*for(int i=0; i< sweeper->get_end_state()->data().size(); i++){
	std::cout << "ergebnis " << sweeper->get_end_state()->data()[i] << " " << sweeper->exact(t_end)->data()[i] <<  std::endl ;
	}*/
	        std::cout <<  "fein" << std::endl;
        auto naeherung = sweeper->get_end_state()->data();
        auto exact     = sweeper->exact(t_end)->data();
        /*for (int i=0; i< sweeper->get_end_state()->data().size(); i++){
          std::cout << sweeper->exact(0)->data()[i] << " " << naeherung[i] << "   " << exact[i] << std::endl;
        }*/

	  typedef Dune::BlockVector<Dune::FieldVector<double, 1> > VectorType;
          VectorType x = sweeper->get_end_state()->data();
	  VectorType y = sweeper->exact(t_end)->data();
        
	  
	  x-=y;
	 // x/=y.infinity_norm();
	  
	  //sweeper->states()[sweeper->get_states().size()-1]->scaled_add(-1.0 , sweeper->exact(t_end));
          //std::cout << sweeper->states()[sweeper->get_states().size()-1]->norm0()<<  std::endl ;

	  sweeper->get_end_state()->scaled_add(-1.0 , sweeper->exact(t_end));
	  std::cout << "FEHLER" << sweeper->get_end_state()->norm0()<< " number solves " << sweeper->num_solves << std::endl ;
	
      ofstream f;
	  stringstream ss;
	  ss << nelements;
	  string s = "solution_sdc/" + ss.str() + ".dat";
	  f.open(s, ios::app | std::ios::out );
      f << nelements << " " << dt << " "<< sweeper->get_end_state()->norm0()<< endl;
	  //f << nelements << " " << dt << " "<< x.infinity_norm()<< endl;

	  f.close();

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

    using sweeper_t = Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM>>;

    pfasst::init(argc, argv, sweeper_t::init_opts);

    const size_t nelements = get_value<size_t>("num_elements", 180); //Anzahl der Elemente pro Dimension
    const size_t nnodes = get_value<size_t>("num_nodes", 3);
    const QuadratureType quad_type = QuadratureType::GaussRadau;
    const double t_0 = 0.0;
    const double dt = get_value<double>("dt", 0.05);
    double t_end = get_value<double>("tend", 0.1);
    size_t nsteps = get_value<size_t>("num_steps", 0);
    double newton = get_value<double>("newton", 0.1);

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
    
    const size_t niter = get_value<size_t>("num_iters", 10);

    pfasst::examples::heat_FE::run_sdc(nelements, BASIS_ORDER, DIM, nnodes, quad_type, t_0, dt, t_end, niter, newton);

  }

#endif
