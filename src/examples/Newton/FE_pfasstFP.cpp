#include <config.h>

#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include <mpi.h>

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>


#include <pfasst/comm/mpi_p2p.hpp>
#include <pfasst/controller/two_level_pfasst.hpp>

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
                      const double& t_0, const double& dt, const double& t_end, const size_t& niter)
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
              fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); //grid->levelGridView(i));//gridView);

        } 


	auto coarse = std::make_shared<SweeperType>(fe_basis[1], 1,  grid);

        auto fine = std::make_shared<SweeperType>(fe_basis[0] , 0, grid);	const auto num_nodes = nnodes;	
    	const auto num_time_steps = 1;//t_end/dt;

	vector<vector<shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t>>>  _new_newton_state_coarse;
	vector<vector<shared_ptr<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>::encap_t>>>  _new_newton_state_fine;	
	//vector<vector<shared_ptr<Dune::BlockVector<Dune::FieldVector<double, 1>>>>>  _new_newton_state_coarse;
	//vector<vector<shared_ptr<Dune::BlockVector<Dune::FieldVector<double, 1>>>>>  _new_newton_state_fine;
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



for(int time=0; time<((t_end-t_0)/dt); time+=num_pro){	
    for(int ne=0; ne<4; ne++){



	TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
        pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);

        

        auto coarse = std::make_shared<SweeperType>(fe_basis[1], 1,  grid);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<SweeperType>(fe_basis[0], 0,  grid);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        
        
        coarse->is_coarse=true;
        fine->is_coarse=false;
                        
        dunetransfer = std::make_shared<TransferOperatorAssembler<GridType>>(*grid);
	transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	for (int i=0; i< n_levels-1; i++){
	      transferMatrix->push_back(new MatrixType()); // hier nur referenz die evtl geloescht wird??
	}
	dunetransfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	    
	std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;

	for (int i=0; i< vecvec->at(0)->N(); i++){
	      for (int j=0; j< (*vecvec->at(0)).M(); j++){
		if(vecvec->at(0)->exists(i,j)){
		}
	      }
        }
        
       
        auto transfer = std::make_shared<TransferType>();
	transfer->create(vecvec);


        fine->set_abs_residual_tol(1e-12);
        coarse->set_abs_residual_tol(1e-12);



	pfasst.add_sweeper(coarse, true);
	pfasst.add_sweeper(fine);
        
        pfasst.add_transfer(transfer);
        std::cout << "nach add ransfer"<< std::endl;

        pfasst.set_options();



        pfasst.status()->time() =  t_0 + time*dt*num_pro;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_0 + (time+1)*dt*num_pro;
        pfasst.status()->max_iterations() = niter;
        std::cout << "******************************** pfasst t0 " << pfasst.status()->time() << "tend " << pfasst.status()->t_end() << std::endl;
        pfasst.setup();

        coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());
        fine->initial_state() = fine->exact(pfasst.get_status()->get_time());


	if(time=0 && ne==0) 	
	for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1; j++){
		for(int k=0; k< _new_newton_state_coarse[i][j]->data().size(); k++){
    		 (*_new_newton_state_coarse[i][j]).data()[k]= 0; 
    		}
		for(int k=0; k< _new_newton_state_fine[i][j]->data().size(); k++){
    		 (*_new_newton_state_fine[i][j]).data()[k]= 0; 
    		}
		}
	}



	for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1; j++){
			//for(int k=0; k< _new_newton_state_coarse[i][j]->data().size(); k++){
    			//coarse->last_newton_state()[i][j]->data()[k] = _new_newton_state_coarse[i][j]->data()[k]  ;
			//}
    			transfer->restrict_u(_new_newton_state_fine[i][j], coarse->last_newton_state()[i][j]);
			for(int k=0; k< _new_newton_state_fine[i][j]->data().size(); k++){
    			fine->last_newton_state()[i][j]->data()[k] = _new_newton_state_fine[i][j]->data()[k]  ;
			}
    		}
	}

	    for(int m=0; m< num_nodes +1; m++){
	    	fine->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(fine->M_dune); 
            	fine->evaluate_df2(*fine->df_dune[0][m], fine->last_newton_state()[0][m]);
	    	coarse->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(coarse->M_dune); 
	    	transfer->restrict_dune_matrix(*fine->df_dune[0][m], *coarse->df_dune[0][m]);
		auto result = fine->get_encap_factory().create();
            	result->zero();
                fine->evaluate_f2(result, fine->last_newton_state()[0][m]);
		fine->df_dune[0][m]->mmv(fine->last_newton_state()[0][m]->data(), result->data());

	    	fine->coarse_rhs()[0][m]->data() =result->data();
		transfer->restrict_data(fine->coarse_rhs()[0][m], coarse->coarse_rhs()[0][m]);
                
	    	//Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > vgl_M = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(coarse->M_dune); ///////M
		//transfer->restrict_dune_matrix(*fine->df_dune[0][m], vgl_M);


	    }




        std::cout << "vor run"<< std::endl;
        pfasst.run();
        pfasst.post_run();

        
        MPI_Barrier(MPI_COMM_WORLD);

        
        if(my_rank==0) {
        auto anfang    = fine->exact(0)->data();
        auto naeherung = fine->get_end_state()->data();
        auto exact     = fine->exact(t_end)->data();
        for (int i=0; i< fine->get_end_state()->data().size(); i++){
          std::cout << anfang[i] << " " << naeherung[i] << "   " << exact[i] << " "  <<  std::endl;
        }

        std::cout << "******************************************* " << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        fine->get_end_state()->scaled_add(-1.0, fine->exact(t_end));
        std::cout << "Fehler: "  << fine->get_end_state()->norm0() << " " << std::endl;
        std::cout << "Time: "  << time << " " << "Newton: " << ne << std::endl;

	std::cout << "******************************************* " << std::endl;
	}

        MPI_Barrier(MPI_COMM_WORLD);
        
                /*if(my_rank==1) {
        auto anfang    = fine->exact(0)->data();
        auto naeherung = fine->get_end_state()->data();
        auto exact     = fine->exact(t_end)->data();
        for (int i=0; i< fine->get_end_state()->data().size(); i++){
          std::cout << anfang[i] << " " << naeherung[i] << "   " << exact[i] << " "  <<  std::endl;
        }

        std::cout << "******************************************* " << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        fine->get_end_state()->scaled_add(-1.0, fine->exact(t_end));
        std::cout << "Fehler: "  << fine->get_end_state()->norm0() << " " << std::endl;


	std::cout << "******************************************* " << std::endl;
	}*/

    	for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1 ; j++){
			for(int k=0; k< _new_newton_state_coarse[i][j]->data().size(); k++)
    				(*_new_newton_state_coarse[i][j]).data()[k] = coarse->new_newton_state()[i][j]->data()[k];
    		
			for(int k=0; k< _new_newton_state_fine[i][j]->data().size(); k++)
    				(*_new_newton_state_fine[i][j]).data()[k] = fine->new_newton_state()[i][j]->data()[k];
    		}
	}
	

        MPI_Barrier(MPI_COMM_WORLD);
}
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


  const size_t nelements = get_value<size_t>("num_elements", 4); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.1);
  double t_end = get_value<double>("tend", 0.4);
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

  run_pfasst(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter);

  pfasst::Status<double>::free_mpi_datatype();

  std::cout << "my rank " << my_rank<< " of "<< num_pro << std::endl;

  MPI_Finalize();

  return 0;
}
