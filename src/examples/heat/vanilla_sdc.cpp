/**
 * Advection-Diffusion with vanilla SDC.
 *
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/vanilla_sdc.cpp
 * @since v0.1.0
 */

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <cstdlib>
#include <memory>

#include <pfasst.hpp>
#include <pfasst/controller/sdc.hpp>
#include "dune_vec.hpp"

#include "assemble.hpp"
#include "advection_diffusion_sweeper.hpp"



namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Advection/diffusion example using an encapsulated IMEX sweeper.
       *
       * This example uses a vanilla SDC sweeper.
       *
       * @ingroup AdvectionDiffusion
       */
      error_map run_vanilla_sdc(double abs_residual_tol, double rel_residual_tol=0.0)
      {
        SDC<> sdc;

        auto const nnodes    = config::get_value<size_t>("num_nodes", 3); // call with ./dune-PFASST num_nodes=4
        auto const nelements = config::get_value<size_t>("num_elem", 64); //number of finite elements per dimension
        auto const quad_type = quadrature::QuadratureType::GaussLegendre;

	
	auto FinEl   = make_shared<fe_manager>(nelements); //spatial_dofs
	auto ndofs   = FinEl->get_ndofs();
        auto quad    = quadrature::quadrature_factory(nnodes, quad_type);
        auto factory = make_shared<encap::Dune_VectorFactory<double>>(ndofs);
        auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(FinEl, ndofs);

        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);
        sweeper->set_residual_tolerances(abs_residual_tol, rel_residual_tol);

        sdc.add_level(sweeper);
        sdc.set_duration(0.0, 1000*0.001, 0.001, 4);
        sdc.setup();
	
	auto q0 = sweeper->get_start_state();	

        sweeper->exact(q0, 0.0);

        sdc.run();

	

	
	typedef encap::Dune_VectorEncapsulation<double, double> vector_type;
        vector_type neu(ndofs); //= encap::as_vector<double,double>(q0);
	
	sweeper->exact(neu, 1.0);	
	//vector_type exact=encap::as_vector<double,double>(neu);
	vector_type u=  encap::as_vector<double,double>(sweeper->get_end_state());
	vector_type u0(ndofs);
	sweeper->exact(u0, 0.0);
	  
	
	//vector_type neu = encap::as_vector<double,double>(q0);	
	
	
	typedef Dune::YaspGrid<1> GridType; 
        auto grid = FinEl->get_grid();
        typedef GridType::LeafGridView GridView;
        GridType::LeafGridView gridView = grid->leafGridView();
        Dune::VTKWriter<GridView> vtkWriter(gridView);
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > VectorType;
	//VectorType x(FinEl->get_ndofs());
	//VectorType x(FinEl->get_ndofs());
	for(int i=0; i< FinEl->get_ndofs(); i++){
	  //std::cout<< u0[i][0] << " " << u[i][0] << " " << neu[i][0]  <<std::endl;
	}
	std::cout << (neu.axpy(-1,u)).infinity_norm() << std::endl; 
        //VectorType y = sweeper->exact(t_end)->data();
        //VectorType z = sweeper->initial_state()->data();
        vtkWriter.addVertexData(neu  , "fe_solution");
        //vtkWriter.addVertexData(y, "exact_solution");
        //vtkWriter.addVertexData(z, "initial_data");
        vtkWriter.write("heat_result");
	
        return sweeper->get_errors();	

      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::init(argc, argv, pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_logs);
  
  pfasst::examples::advection_diffusion::run_vanilla_sdc(0.0);
}
#endif