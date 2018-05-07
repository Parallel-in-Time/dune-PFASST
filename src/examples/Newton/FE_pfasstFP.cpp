#include <config.h>

#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include <mpi.h>

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>


#include <pfasst/comm/mpi_p2p.hpp>
#include <pfasst/controller/two_level_pfasst_n.hpp>

#include "fischer_sweeper.hpp"
#include "../../datatypes/dune_vec.hpp"
#include "spectral_transfer.hpp"

#include <vector>


using namespace pfasst::examples::fischer_example;


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;
using pfasst::encap::DuneEncapsulation;

using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::contrib::SpectralTransfer;
using pfasst::TwoLevelPfasst;
typedef pfasst::comm::MpiP2P CommunicatorType;


typedef DuneEncapsulation<double, double, 1>                     EncapType;


using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  
using SweeperType = fischer_sweeper<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>,   FE_function >;


typedef pfasst::transfer_traits<SweeperType, SweeperType, 1>       TransferTraits;
typedef SpectralTransfer<TransferTraits>                           TransferType;



      void run_pfasst(const size_t nelements, const size_t basisorder, const size_t dim, const size_t& nnodes, const pfasst::quadrature::QuadratureType& quad_type,
                      const double& t_0, const double& dt, const double& t_end, const size_t& niter, double newton)
      {


        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        

                 
        typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
        typedef GridType::LevelGridView GridView;
        using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
    
        std::shared_ptr<TransferOperatorAssembler<GridType>> dunetransfer;

        std::shared_ptr<std::vector<MatrixType*>> transferMatrix;

        std::shared_ptr<GridType> grid;

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
              fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); 

        } 


	auto coarse = std::make_shared<SweeperType>(fe_basis[1], 1,  grid);


        auto fine = std::make_shared<SweeperType>(fe_basis[0] , 0, grid);	const auto num_nodes = nnodes;	
    	const auto num_time_steps = 1;//t_end/dt;

	vector<vector<shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t>>>  _new_newton_state_coarse;
	vector<vector<shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t>>>  _new_newton_state_fine;

	shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t> _new_initial_coarse;
	shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t> _new_initial_fine;
	shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t> _copy_end_state;

	_copy_end_state = fine->get_encap_factory().create();
	_new_initial_coarse = coarse->get_encap_factory().create();
	_new_initial_fine = fine->get_encap_factory().create();
	

    	_new_newton_state_coarse.resize(num_time_steps);
    	_new_newton_state_fine.resize(num_time_steps);

    	for(int i=0; i< num_time_steps; i++){	
		_new_newton_state_fine[i].resize(num_nodes + 1);
		_new_newton_state_coarse[i].resize(num_nodes + 1);
		for(int j=0; j<num_nodes +1 ; j++){
			_new_newton_state_fine[i][j] =  fine->get_encap_factory().create(); //std::make_shared<Dune::BlockVector<Dune::FieldVector<double, 1>>>(fe_basis[0]->size());
			_new_newton_state_coarse[i][j] = coarse->get_encap_factory().create(); //std::make_shared<Dune::BlockVector<Dune::FieldVector<double, 1>>>(fe_basis[1]->size());
		}
    	}




	/*auto coarse_initial = std::make_shared<SweeperType>(fe_basis[1], 1,  grid); 

	coarse_initial->quadrature() = quadrature_factory<double>(nnodes, quad_type);

	auto sweeper = std::make_shared<sweeper_t>(fe_basis[1] , 1, grid); 
    	auto sdc = std::make_shared<heat_FE_sdc_t>();    


	sweeper->is_coarse = false;
    	sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);

    	sdc->add_sweeper(sweeper);
    	sdc->set_options();
    	if (my_rank==0){ sdc->status()->time() = 0;}else{sdc->status()->time() = 0.5;}
    	sdc->status()->dt() = 0.5;
    	if (my_rank==0){sdc->status()->t_end() = 0.5;}else{sdc->status()->t_end() = 1;}
    	sdc->status()->max_iterations() = 2;
    	sdc->setup();

        if (my_rank==0){ 

	sweeper->initial_state() = sweeper->exact(0);
	
        sdc->run();

        sdc->post_run();
	}
	MPI_Barrier(MPI_COMM_WORLD);*/



std::cout << "num_pro " << num_pro << std::endl;
int num_solves=0;
for(int time=0; time<((t_end-t_0)/dt); time+=num_pro){	

    std::cout << "-------------------------------------------------------------------------------------------------------  time " << time << " " << std::endl;	
    for(int ne=0; ne<12; ne++){

	TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
        pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);

        auto coarse = std::make_shared<SweeperType>(fe_basis[1], 1,  grid);
	auto coarse_helper = std::make_shared<SweeperType>(fe_basis[1], 1,  grid); 
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        coarse_helper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<SweeperType>(fe_basis[0], 0,  grid);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        coarse->is_coarse=true;
	coarse_helper->is_coarse=true;
        fine->is_coarse=false;
                        
        dunetransfer = std::make_shared<TransferOperatorAssembler<GridType>>(*grid);
	transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	for (int i=0; i< n_levels-1; i++){
	      transferMatrix->push_back(new MatrixType()); 
	}
	dunetransfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	    
	std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;

        
       
        auto transfer = std::make_shared<TransferType>();
	transfer->create(vecvec);

	fine->num_solves+=num_solves;

	

	pfasst.add_sweeper(coarse, true);
	pfasst.add_sweeper(fine);
        
        pfasst.add_transfer(transfer);
        pfasst.set_options();

        pfasst.status()->time() =  t_0 + time*dt;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_0 + (time+num_pro)*dt;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

	//if(time==0){
   
	//}else {
	//	for (int i=0; i< fine->initial_state()->data().size(); i++) fine->initial_state()->data()[i] = _new_initial_fine->data()[i];   
	//}


	if(time==0 && ne==0){ 	
		//for(int i=0; i< num_time_steps; i++){	
			for(int j=0; j<num_nodes +1; j++){
				/*for(int k=0; k<(_new_newton_state_coarse[i][j])->data().size(); k++){
					(_new_newton_state_coarse[i][j])->data()[k] = fine->exact( pfasst.status()->time()); //0;
				}
				for(int k=0; k<(_new_newton_state_fine[i][j])->get_data().size(); k++){
					(_new_newton_state_fine[i][j])->data()[k] = coarse->exact( pfasst.status()->time());
				}*/
				(_new_newton_state_fine[0][j])  = fine->exact( pfasst.status()->time());
				(_new_newton_state_coarse[0][j]) = coarse->exact( pfasst.status()->time());
				for(int k=0; k<(_new_newton_state_fine[0][j])->data().size(); k++){
					if(my_rank==0) std::cout << "start fine" <<  (_new_newton_state_fine[0][j])->data()[k] << std::endl;
				}
			}

		//}
	}


        fine->initial_state() = (_new_newton_state_fine[0][0]);     
        coarse->initial_state() = (_new_newton_state_coarse[0][0]);  
	//MPI_Bcast(&fine->initial_state()->data()[0], fine->initial_state()->data().size(), MPI_FLOAT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&coarse->initial_state()->data()[0], coarse->initial_state()->data().size(), MPI_FLOAT, 0, MPI_COMM_WORLD);

	/*MPI_Status Stat;
	if(my_rank==0) 	MPI_Send(&(_new_newton_state_fine[0][_new_newton_state_fine[0][0]->data().size()-1]), _new_newton_state_fine[0][0]->data().size(),  MPI_FLOAT, 1,7, MPI_COMM_WORLD);
	if(my_rank==1) 	MPI_Recv(&(fine->initial_state()->data()[0]), _new_newton_state_fine[0][0]->data().size(),  MPI_FLOAT, 0,7, MPI_COMM_WORLD, &Stat);

	if(my_rank==0) 	MPI_Send(&(_new_newton_state_fine[0][_new_newton_state_coarse[0][0]->data().size()-1]), _new_newton_state_coarse[0][0]->data().size(),  MPI_FLOAT, 1,7, MPI_COMM_WORLD);
	if(my_rank==1) 	MPI_Recv(&(coarse->initial_state()->data()[0]), _new_newton_state_coarse[0][0]->data().size(),  MPI_FLOAT, 0,7, MPI_COMM_WORLD, &Stat);*/

 	/*for (size_t predict_step = 0;
         predict_step <= this->get_communicator()->get_rank();
         ++predict_step) {
      	// do the sweeper's prediction once ...
      	if (my_rank ) {
        //this->predict_coarse();
	//this->get_coarse()->spread(0);
      	} else {
        // and default sweeps for subsequent processes
        //this->recv_coarse();
        //this->sweep_coarse();
      	}*/


	//for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1; j++){
    			transfer->restrict_u(_new_newton_state_fine[0][j], coarse->last_newton_state()[0][j]);
			for(int k=0; k< _new_newton_state_fine[0][j]->data().size(); k++){
    				fine->last_newton_state()[0][j]->data()[k] = _new_newton_state_fine[0][j]->data()[k]  ;
    				fine->new_newton_state()[0][j]->data()[k] = _new_newton_state_fine[0][j]->data()[k]  ;
    				//coarse->last_newton_state()[i][j]->data()[k] = _new_newton_state_coarse[i][j]->data()[k]  ;
			}
			for(int k=0; k< _new_newton_state_coarse[0][j]->data().size(); k++){
    				//fine->last_newton_state()[i][j]->data()[k] = _new_newton_state_fine[i][j]->data()[k]  ;
    				coarse->last_newton_state()[0][j]->data()[k] = _new_newton_state_coarse[0][j]->data()[k]  ;
    				coarse->new_newton_state()[0][j]->data()[k] = _new_newton_state_coarse[0][j]->data()[k]  ;
			}
			//if (my_rank==1) fine->last_newton_state()[i][j] =  fine->exact(1) ;
			//if (my_rank==0) fine->last_newton_state()[i][j] =  fine->exact(0.5) ;
    		}
	//}












	    for(int m=0; m< num_nodes +1; m++){
	    	fine->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(fine->M_dune); 
            	fine->evaluate_df2(*fine->df_dune[0][m], fine->last_newton_state()[0][m]);

	    	coarse->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(coarse->M_dune); 
	    	transfer->restrict_dune_matrix(*fine->df_dune[0][m], *coarse->df_dune[0][m]);
		
		auto result = fine->get_encap_factory().create();
            	result->zero();
                fine->evaluate_f2(result, fine->last_newton_state()[0][m]);
		fine->df_dune[0][m]->mmv(fine->last_newton_state()[0][m]->data(), result->data());
		//fine->evaluate_f3(result, fine->last_newton_state()[0][m]);	    	
		//auto mat = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(fine->M_dune); 
		//fine->evaluate_df3(*mat, fine->last_newton_state()[0][m]);
		//mat->mv(fine->last_newton_state()[0][m], result->data());

		fine->coarse_rhs()[0][m]->data() =result->data();
		transfer->restrict_data(fine->coarse_rhs()[0][m], coarse->coarse_rhs()[0][m]);
	    }


        //if(ne ==3) fine->set_abs_residual_tol(1e-1);
        //if(ne ==3) coarse->set_abs_residual_tol(1e-1);

        pfasst.run();
        pfasst.post_run();

        
	std::cout << my_rank << " nach dem run " << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        
        if(my_rank==num_pro-1) {
        auto anfang    = fine->exact(0)->data();
        auto naeherung = fine->get_end_state()->data();
        auto exact     = fine->exact( t_0 + (time+num_pro)*dt)->data();
        auto exact2     = fine->exact( t_0 + (time+1)*dt)->data();
        for (int i=0; i< fine->get_end_state()->data().size(); i++){
          std::cout <<  t_0 + (time+num_pro)*dt << anfang[i] << " " << naeherung[i] << "   " << exact[i] << " "  <<  std::endl;
        }
	}





	

	for (int i=0; i< fine->get_end_state()->data().size(); i++) _copy_end_state->data()[i] = fine->get_end_state()->data()[i];
	fine->get_end_state()->scaled_add(-1.0, _new_newton_state_fine[num_time_steps-1][num_nodes]);
        std::cout << my_rank << " NEWTON *****************************************      Fehler: "  << fine->get_end_state()->norm0() << " " << std::endl;
	std::cout << my_rank << " num_solves  " << fine->num_solves <<  std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

	int local=0, global=0;
	if (fine->get_end_state()->norm0()<newton) local=1;	
	MPI_Reduce(&local, &global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&global, 1, MPI_INT, 0, MPI_COMM_WORLD);



	if(global>0 || ne==3) {
		//if (my_rank == num_pro-1){for (int i=0; i< fine->get_end_state()->data().size(); i++) _new_initial_fine->data()[i] = _copy_end_state->data()[i];}
		//MPI_Bcast(&(_new_initial_fine->data()[0]),_new_initial_fine->data().size(), MPI_DOUBLE, num_pro-1,MPI_COMM_WORLD);



		/*if(t_0 + (time+num_pro)*dt == t_end){	
			int final_solves=0;				
			MPI_Reduce(&fine->num_solves, &final_solves, 1, MPI_INT, MPI_SUM, num_pro-1, MPI_COMM_WORLD);
        		if (my_rank == num_pro-1){ std::cout << "NEWTON *****************************************      Fehler: "  << fine->get_end_state()->norm0() << " " << std::endl;         	std::cout << "Fehler am Ender : "  << _copy_end_state->norm0() << " " << std::endl;

        			ofstream f;
        			stringstream ss;
        			ss << nelements;
        			string s = "solution_pfasst12/" + ss.str() + ".dat";
        			f.open(s, ios::app | std::ios::out );
        			f << nelements << " " << dt << " "<< _copy_end_state->norm0() << " number solves " << final_solves << endl;
        			f.close();
        			std::cout << "******************************************* " << std::endl;
			}
        	


		}*/
	fine->new_newton_state()[0][num_nodes]->scaled_add(-1.0, fine->exact( t_0 + (1+my_rank)*dt));
       	std::cout << my_rank << " Fehler am Ender : "  << fine->new_newton_state()[0][num_nodes]->norm0() << " " << std::endl;

		std::cout << " vor dem break " << std::endl; break;
	}



	num_solves = fine->num_solves;
        MPI_Barrier(MPI_COMM_WORLD);

	std::cout << "rank " << my_rank<< std::endl;

    	//for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1 ; j++){
			for(int k=0; k< _new_newton_state_coarse[0][j]->data().size(); k++)
    				(_new_newton_state_coarse[0][j])->data()[k] = coarse->states()[j]->data()[k]; //coarse->new_newton_state()[0][j]->data()[k];
    		
			for(int k=0; k< _new_newton_state_fine[0][j]->data().size(); k++)
    				(_new_newton_state_fine[0][j])->data()[k] = fine->states()[j]->data()[k]; //new_newton_state()[0][j]->data()[k];
    		}
	//}
	

	fine->new_newton_state()[0][num_nodes]->scaled_add(-1.0, fine->exact( t_0 + (1+my_rank)*dt));
       	std::cout << my_rank << " Fehler am Ender : "  << fine->new_newton_state()[0][num_nodes]->norm0() << " " << std::endl;


        MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "after last barrier rank " << my_rank<< std::endl;


}

std::cout << "-------------------------------------------------------------------------------------------------------  beginn weiterer zeitschritt" << time << " " << std::endl;
}


      }



int main(int argc, char** argv)
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  MPI_Init(&argc, &argv);
        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

        
  
  
  pfasst::init(argc, argv, SweeperType::init_opts);
  pfasst::Status<double>::create_mpi_datatype();


  const size_t nelements = get_value<size_t>("num_elements", 32); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0;
  const double dt = get_value<double>("dt", 0.1);
  double t_end = get_value<double>("tend", 0.2);
  size_t nsteps = get_value<size_t>("num_steps", 0);
  double newton = get_value<double>("newton", 0);
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
  const size_t niter = get_value<size_t>("num_iters", 50);

  run_pfasst(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter, newton);

  pfasst::Status<double>::free_mpi_datatype();

  MPI_Barrier(MPI_COMM_WORLD);



  MPI_Finalize();

  return 0;
}
