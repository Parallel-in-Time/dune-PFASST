#include <config.h>

#include <fenv.h>

#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/two_level_mlsdc.hpp>




#include <vector>


//#include "FE_sweeper.hpp"
#include "fischer_sweeper.hpp"
#include "../../datatypes/dune_vec.hpp"
#include "spectral_transfer.hpp"


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
                                           const size_t niter) {



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

        auto fine = std::make_shared<sweeper_t_coarse>(fe_basis[0] , 0, grid);


	const auto num_nodes = nnodes;	
    	const auto num_time_steps = t_end/dt;

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




    for(int ne=0; ne<4; ne++){

        auto mlsdc = std::make_shared<heat_FE_mlsdc_t>();

        auto coarse = std::make_shared<sweeper_t_coarse>(fe_basis[1], 1,  grid);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);


        auto fine = std::make_shared<sweeper_t_coarse>(fe_basis[0] , 0, grid);
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
	    //std::cout <<  "transfer erzeugt groesse " << (*vecvec->at(0)).M() <<  std::endl;
	for (int i=0; i< vecvec->at(0)->N(); i++){
	      for (int j=0; j< (*vecvec->at(0)).M(); j++){
		if(vecvec->at(0)->exists(i,j)){
		  //std::cout << ((*vecvec->at(0))[i][j]) << std::endl;
		}
	      }
        }
        
        auto transfer = std::make_shared<transfer_t>();
        transfer->create(vecvec);
           
        
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



	

	if(ne==0) 	
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


	Dune::BlockVector<Dune::FieldVector<double,1>> rM_rv;
	Dune::BlockVector<Dune::FieldVector<double,1>> rMv;
	Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > vgl_M = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(coarse->M_dune); ///////M
	transfer->restrict_dune_matrix(fine->M_dune, vgl_M);
	
        /*std::cout << "coarse Stiffnessmatrix" << std::endl;
	for(int i=0; i< coarse->M_dune.M(); i++)
		for(int j=0; j< coarse->M_dune.M(); j++)
			if(vgl_M.exists(i,j)) std::cout << "coarse " << coarse->M_dune[i][j] << "restringiertes coarse "<< vgl_M[i][j] << std::endl;
	if (ne==0) std::exit(0);*/
           



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


	
		



            









        std::cout << "starting run" << std::endl;

        
        


        mlsdc->run();
        
        mlsdc->post_run();



        std::cout << "******************************************* " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << "Fehler: " <<  std::endl ;
        fine->states()[fine->get_states().size()-1]->scaled_add(-1.0 , fine->exact(t_end));
        std::cout << fine->states()[fine->get_states().size()-1]->norm0()<<  std::endl ;
        std::cout << "number states " << fine->get_states().size() << std::endl ;
        std::cout << "******************************************* " <<  std::endl ;


    	for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1 ; j++){
		for(int k=0; k< _new_newton_state_coarse[i][j]->data().size(); k++)
    		(*_new_newton_state_coarse[i][j]).data()[k] = coarse->new_newton_state()[i][j]->data()[k];
    		
		for(int k=0; k< _new_newton_state_fine[i][j]->data().size(); k++)
    		(*_new_newton_state_fine[i][j]).data()[k] = fine->new_newton_state()[i][j]->data()[k];
    		}
		}
	}



       

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

  const size_t nelements = get_value<size_t>("num_elements", 180); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  //const size_t ndofs = get_value<size_t>("num_dofs", 8);
  const size_t coarse_factor = get_value<size_t>("coarse_factor", 1);
  //const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.05);
  double t_end = get_value<double>("tend", 0.05);
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

  run_mlsdc(nelements, BASE_ORDER, DIMENSION, coarse_factor, nnodes, quad_type, t_0, dt, t_end, niter);
}
#endif 
