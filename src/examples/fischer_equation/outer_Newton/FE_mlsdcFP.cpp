#include <config.h>

#include <fenv.h>

#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/two_level_mlsdc_n.hpp>




#include <vector>


#include "../1d_transfer/fe_manager.hpp"
#include "fischer_sweeper.hpp"
#include <pfasst/encap/dune_vec.hpp>
#include "../1d_transfer/spectral_transfer.hpp"


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;

using namespace pfasst::examples::fischer_example;




      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::TwoLevelMLSDC;
      using pfasst::quadrature::QuadratureType;

      
      using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  
      using sweeper_t_coarse = fischer_sweeper<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>,   FE_function >;
      using transfer_traits_t = pfasst::transfer_traits<sweeper_t_coarse, sweeper_t_coarse, 1>;
      using transfer_t = SpectralTransfer<transfer_traits_t>;
      using heat_FE_mlsdc_t = TwoLevelMLSDC<transfer_t>;


      void run_mlsdc(const size_t nelements, const size_t basisorder, const size_t DIM, const size_t coarse_factor,
                                           const size_t nnodes, const QuadratureType& quad_type,
                                           const double& t_0, const double& dt, const double& t_end,
                                           const size_t niter, double newton, double tol) {



//         auto FinEl = make_shared<fe_manager>(nelements, 2);

        
                
        typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
        typedef GridType::LevelGridView GridView;
        using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
    
        std::shared_ptr<TransferOperatorAssembler<GridType>> dunetransfer;

        std::shared_ptr<std::vector<MatrixType*>> transferMatrix;

        std::shared_ptr<GridType> grid;

        int n_levels=2;

        std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); ; 
        //std::vector<std::shared_ptr<BasisFunction> > fe_basis_p;

    
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
	      //n_dof[n_levels-i-1]    = fe_basis[n_levels-i-1]->size();
        } 
        
        
        
        using pfasst::quadrature::quadrature_factory;


	auto coarse = std::make_shared<sweeper_t_coarse>(fe_basis[1], 1,  grid);

        auto fine = std::make_shared<sweeper_t_coarse>(fe_basis[0] , 0, grid);    //[0]


    
	const auto num_nodes = nnodes;	
    	const auto num_time_steps = 1; //t_end/dt;

	vector<shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t>>  _new_newton_state_coarse;
	vector<shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t>>  _new_newton_state_fine;	



	_new_newton_state_fine.resize(num_nodes + 1);
	_new_newton_state_coarse.resize(num_nodes + 1);
	for(int j=0; j<num_nodes +1 ; j++){
			_new_newton_state_fine[j] =  fine->get_encap_factory().create(); //std::make_shared<Dune::BlockVector<Dune::FieldVector<double, 1>>>(fe_basis[0]->size());
			_new_newton_state_coarse[j] = coarse->get_encap_factory().create(); //std::make_shared<Dune::BlockVector<Dune::FieldVector<double, 1>>>(fe_basis[1]->size());
	}
   
	Dune::BlockVector<Dune::FieldVector<double, 1>> _new_newton_initial_coarse(fe_basis[1]->size());    
	Dune::BlockVector<Dune::FieldVector<double, 1>> _new_newton_initial_fine(fe_basis[0]->size()); 
    	
	std::cout.precision ( 10 );

for(int time=0; time<(t_end-t_0)/dt; time++){	

    for(int ne=0; ne<10; ne++){

        auto mlsdc = std::make_shared<heat_FE_mlsdc_t>();

        auto coarse = std::make_shared<sweeper_t_coarse>(fe_basis[1], 1,  grid);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);


        auto fine = std::make_shared<sweeper_t_coarse>(fe_basis[0] , 0, grid);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);


        coarse->is_coarse=true;
        fine->is_coarse=false;
        
        
        fine->set_abs_residual_tol(tol);
        coarse->set_abs_residual_tol(tol);

        std::cout << "erstelle tranfer" << std::endl;
        
        dunetransfer = std::make_shared<TransferOperatorAssembler<GridType>>(*grid);
	transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	for (int i=0; i< n_levels-1; i++){
	      transferMatrix->push_back(new MatrixType()); 
	}
	dunetransfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	    
	std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;


        
        auto transfer = std::make_shared<transfer_t>();
        transfer->create(vecvec);
           
        std::cout << "nach erstelle tranfer" << std::endl;
        
        mlsdc->add_sweeper(coarse, true);
	mlsdc->add_sweeper(fine);

        mlsdc->add_transfer(transfer);


        mlsdc->set_options();



        mlsdc->status()->time() = t_0 + time*dt;
        mlsdc->status()->dt() = dt;
        mlsdc->status()->t_end() = t_0 + (time+1)*dt; 
        mlsdc->status()->max_iterations() = niter;



        mlsdc->setup();
        
        if (time==0){
        	coarse->initial_state() = coarse->exact(mlsdc->get_status()->get_time());
        	fine->initial_state() = fine->exact(mlsdc->get_status()->get_time());
	}else{
		coarse->initial_state()->data() = _new_newton_initial_coarse; 
		fine->initial_state()->data() = _new_newton_initial_fine; 
	}




	

	if(time==0 && ne==0) 	
		for(int j=0; j<num_nodes +1; j++){
    		 	(*_new_newton_state_coarse[j]) = coarse->exact(mlsdc->get_status()->get_time())->data(); 
    		 	(*_new_newton_state_fine[j]) = fine->exact(mlsdc->get_status()->get_time())->data(); 

		}



	for(int j=0; j<num_nodes +1; j++){
			for(int k=0; k< _new_newton_state_coarse[j]->data().size(); k++){
    				coarse->last_newton_state()[j]->data()[k] = _new_newton_state_coarse[j]->data()[k]  ;
			}


			for(int k=0; k< _new_newton_state_fine[j]->data().size(); k++){
    				fine->last_newton_state()[j]->data()[k] = _new_newton_state_fine[j]->data()[k]  ;
			}
    	}
	



           


	    for(int m=0; m< num_nodes +1; m++){
	    	fine->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(fine->M_dune); 
            	fine->evaluate_df2(*fine->df_dune[0][m], fine->last_newton_state()[m]);

	    	coarse->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(coarse->M_dune); 
            	coarse->evaluate_df2(*coarse->df_dune[0][m], coarse->last_newton_state()[m]);	    	


		auto result = fine->get_encap_factory().create();
            	result->zero();

                fine->evaluate_f2(result, fine->last_newton_state()[m]);
		fine->df_dune[0][m]->mmv(fine->last_newton_state()[m]->data(), result->data());

	    	fine->coarse_rhs()[0][m]->data() =result->data();
		

		auto resultc = coarse->get_encap_factory().create();
            	resultc->zero();

                coarse->evaluate_f2(resultc, coarse->last_newton_state()[m]);
		coarse->df_dune[0][m]->mmv(coarse->last_newton_state()[m]->data(), resultc->data());

	    	coarse->coarse_rhs()[0][m]->data() =resultc->data();



	    }




        std::cout << "starting run" << std::endl;

        
        


        mlsdc->run();
        
        mlsdc->post_run();

        auto anfang    = fine->exact(0)->data();
        auto naeherung = fine->get_end_state()->data();
        auto exact     = fine->exact( t_0 + (time+1)*dt)->data();


        std::cout << "******************************************* " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << "Fehler zur Zeit  " << t_0 + (time+1)*dt << std::endl ;
        fine->get_end_state()->scaled_add(-1.0 , fine->exact(t_0 + (time+1)*dt));
        std::cout << fine->get_end_state()->norm0()<<  std::endl ;
        std::cout << "time step " << time << std::endl ;
        std::cout << "******************************************* " <<  std::endl ;
	std::cout << "groesse loesungsvektor " << fine->get_end_state()->data().size() << std::endl ;
	std::cout << "Parameter " << fine->_n << " " << fine->_nu << std::endl ;

	(*_new_newton_state_fine[num_nodes]).data() -= fine->new_newton_state()[num_nodes]->data();
	
	//(*_new_newton_state_fine[num_nodes]).data() -= fine->get_end_state()->data();
	
	
	
	
	std::cout << "******************************* Newton " << (_new_newton_state_fine[num_nodes])->norm0() <<  std::endl ;
	if((_new_newton_state_fine[num_nodes])->norm0()< newton){

		for(int j=0; j<num_nodes +1 ; j++){
		for(int k=0; k< _new_newton_state_coarse[j]->data().size(); k++){
    			(*_new_newton_state_coarse[j]).data()[k] = coarse->new_newton_state()[j]->data()[k];
			//std::cout << "coarse newton solution " << coarse->new_newton_state()[i][j]->data()[k] << std::endl;
    		}
		for(int k=0; k< _new_newton_state_fine[j]->data().size(); k++){
    			(*_new_newton_state_fine[j]).data()[k] = fine->new_newton_state()[j]->data()[k];
			//std::cout << "fine newton solution " << fine->new_newton_state()[i][j]->data()[k] << std::endl;//
		}
    		}
		
		
		for(int j=0; j<num_nodes +1 ; j++){
		for(int k=0; k< _new_newton_state_coarse[j]->data().size(); k++){
    			_new_newton_initial_coarse[k] = coarse->new_newton_state()[j]->data()[k];
			//std::cout << "coarse newton solution " << coarse->new_newton_state()[i][j]->data()[k] << std::endl;
    		}
		for(int k=0; k< _new_newton_state_fine[j]->data().size(); k++){
    			_new_newton_initial_fine[k] = fine->new_newton_state()[j]->data()[k];
			//std::cout << "fine newton solution " << fine->new_newton_state()[i][j]->data()[k] << std::endl;//
		}
    		}
		
		
		std::cout << "************************************* STARTING NEW TIMESTEP "<< time << std::endl;
	
		break;}


		for(int j=0; j<num_nodes +1 ; j++){
		for(int k=0; k< _new_newton_state_coarse[j]->data().size(); k++){
    			(*_new_newton_state_coarse[j]).data()[k] = coarse->new_newton_state()[j]->data()[k];
			//std::cout << "coarse newton solution " << coarse->new_newton_state()[i][j]->data()[k] << std::endl;
    		}
		for(int k=0; k< _new_newton_state_fine[j]->data().size(); k++){
    			(*_new_newton_state_fine[j]).data()[k] = fine->new_newton_state()[j]->data()[k];
			//std::cout << "fine newton solution " << fine->new_newton_state()[i][j]->data()[k] << std::endl;//
		}
    		}
	
	
	}//newton




	
	
	}//time



       

      

}


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
   
      Dune::MPIHelper::instance(argc, argv);  
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;
  //using sweeper_t      = pfasst::examples::heat_FE::Heat_FE<pfasst::sweeper_traits<encap_traits_t,1,1>>;
  //using sweeper_t_fine = pfasst::examples::heat_FE::Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, 2, DIMENSION>>;
  using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  
  using sweeper_t_fine = fischer_sweeper<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>,   FE_function >;
  
  pfasst::init(argc, argv, sweeper_t_fine::init_opts);

  const size_t nelements = get_value<size_t>("num_elements", 32); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 2);
  //const size_t ndofs = get_value<size_t>("num_dofs", 8);
  const size_t coarse_factor = get_value<size_t>("coarse_factor", 1);
  //const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.1);
  double t_end = get_value<double>("tend", 0.1);
  size_t nsteps = get_value<size_t>("num_steps", 0);
  const double newton = get_value<double>("newton", 0.1);  
  const double tol = get_value<double>("abs_res_tol", 1e-12);                    // size of timesteping
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

  run_mlsdc(nelements, BASE_ORDER, DIMENSION, coarse_factor, nnodes, quad_type, t_0, dt, t_end, niter, newton, tol);
}
#endif 
