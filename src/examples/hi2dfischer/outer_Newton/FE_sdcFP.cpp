#include <config.h>
#include "../2d_transfer/dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc_n.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include <pfasst/encap/dune_vec.hpp>
#include "fischer_sweeper.hpp"



#include <dune/fufem/functionspacebases/p1nodalbasis.hh>
#include <dune/fufem/assemblers/operatorassembler.hh>
#include <dune/fufem/assemblers/functionalassembler.hh>


using namespace pfasst::examples::fischer_example;

using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>; 


//using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  




using sweeper_t = fischer_sweeper<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>,   BasisFunction >;

using pfasst::transfer_traits;
using pfasst::contrib::SpectralTransfer;
using pfasst::SDC;
using pfasst::quadrature::QuadratureType;
using heat_FE_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

using pfasst::config::get_value;
using pfasst::quadrature::QuadratureType;
using pfasst::quadrature::quadrature_factory;

using pfasst::examples::fischer_example::fischer_sweeper;





int main(int argc, char** argv) {
      
    	Dune::MPIHelper::instance(argc, argv);
    
    	pfasst::init(argc, argv, sweeper_t::init_opts);

    	const size_t nelements = get_value<size_t>("num_elements", 16);     // spacial dimension: number of grid points per dimension on the coase level    
    	const double t_0 = 0.0;                                             // left point of the time intervall is zero 
    	const double dt = get_value<double>("dt", 0.1);                     // size of timesteping
    	double t_end = get_value<double>("tend", 0.1);                      // right point of the time intervall  
    	const size_t nnodes = get_value<size_t>("num_nodes", 3);            // time intervall: number of sdc quadrature points
    	const QuadratureType quad_type = QuadratureType::GaussRadau;        // quadrature type
    	const size_t niter = get_value<size_t>("num_iters", 10);            // maximal number of sdc iterations
    	const double newton = get_value<double>("newton", 0.1);             // size of timesteping
    	bool output = get_value<double>("output", 0);                       // size of timesteping
    

	//start solving procedure und stop time
	MPI_Barrier(MPI_COMM_WORLD);
    	auto st = MPI_Wtime();
    	
    	/*std::shared_ptr<TransferOperatorAssembler<GridType>> transfer;
    	std::shared_ptr<std::vector<MatrixType*>> transferMatrix;
    	std::shared_ptr<GridType> grid;
    	int n_levels=2;
    	std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); ; 
        Dune::FieldVector<typename GridType::ctype,DIMENSION> L;
        L[0]=1; L[1]=1;
        typename std::array<int,DIMENSION> s;
        std::fill(s.begin(), s.end(), nelements);
        std::bitset<DIMENSION> periodic;
        periodic[0]=true;  
        periodic[1]=true; 
        
        //make sure that the grid will not be automaticly splitted
#if HAVE_MPI
        grid        = std::make_shared<GridType>(L,s,periodic,0, MPI_COMM_SELF);	
#else          
        grid        = std::make_shared<GridType>(L,s,periodic,0);	      
#endif


	fe_basis[1] = std::make_shared<BasisFunction>(grid->levelGridView(0));
	
	grid->globalRefine(1);	
	
	fe_basis[0] = std::make_shared<BasisFunction>(grid->levelGridView(1));

    	//for (int i=0; i<n_levels; i++){	      
	//      grid->globalRefine((bool) i);
	//      auto view = grid->levelGridView(i);
        //      fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); 
    	//} */

	auto FinEl   = make_shared<fe_manager>(nelements,2);

	const auto num_nodes = nnodes;	

	
	vector<shared_ptr<Dune::BlockVector<Dune::FieldVector<double, 1>>>>  _new_newton_state; 

	_new_newton_state.resize(num_nodes + 1);
	
	
	for(int j=0; j<num_nodes +1 ; j++){
		_new_newton_state[j] = std::make_shared<Dune::BlockVector<Dune::FieldVector<double, 1>>>(FinEl->get_basis2()->size());
	}

    	
    	Dune::BlockVector<Dune::FieldVector<double, 1>> _new_initial_state(FinEl->get_basis2()->size());

	//Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > A;
	//Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M;
	//assembleProblem(fe_basis[0], A, M);

	int num_solves = 0;
	
	for(int time=0; time<(t_end-t_0)/dt; time++){	//Zeitschritte

		
    		for(int ne=0; ne<10; ne++){	//Newtonschritte


			auto sweeper = std::make_shared<sweeper_t>(FinEl->get_basis2() , 0, FinEl->get_grid()); 
			sweeper->is_coarse = false;
    			sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
    			auto sdc = std::make_shared<heat_FE_sdc_t>();    
    			sdc->add_sweeper(sweeper);
    			sdc->set_options();
    			sdc->status()->time() = t_0 + time*dt;
    			sdc->status()->dt() = dt;
    			sdc->status()->t_end() = t_0 + (time+1)*dt;
    			sdc->status()->max_iterations() = niter;
    			sdc->setup();
			sweeper->num_solves+=num_solves;
			

			if(time==0 ) {	//im ersten Newton Lauf Anfangswerte setzen
			
        			sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());			
			
				if(ne==0)
					for(int j=0; j<num_nodes +1; j++){
    		 				(*_new_newton_state[j]) = sweeper->initial_state()->data() ;
					}
			


			}else{
				sweeper->initial_state()->data() = _new_initial_state; 
			}


			for(int j=0; j<num_nodes +1; j++){
				for(int k=0; k< _new_newton_state[j]->size(); k++){
    					(sweeper->last_newton_state()[j])->data()[k] = (*_new_newton_state[j])[k]  ;
				}
    			}



			for(int m=0; m< num_nodes +1; m++){
	    			sweeper->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(sweeper->M_dune); 
            			sweeper->evaluate_df2(*sweeper->df_dune[0][m], sweeper->last_newton_state()[m]);
				auto result = sweeper->get_encap_factory().create();
            			result->zero();
                		sweeper->evaluate_f2(result, sweeper->last_newton_state()[m]);
				sweeper->df_dune[0][m]->mmv(sweeper->last_newton_state()[m]->data(), result->data());
	
	    			sweeper->coarse_rhs()[0][m]->data() =result->data();  


			}


    			sdc->run();   
			sdc->post_run();
			num_solves = sweeper->num_solves;


			//does Newton already converge? check residuum:
			(*_new_newton_state[num_nodes]) -= sweeper->get_end_state()->data();
			
			
        		std::cout << "NEWTON *****************************************      Residuum: "  << (*_new_newton_state[num_nodes]).infinity_norm() << " " << std::endl;
			std::cout << "################################################################  num_solves  " << sweeper->num_solves <<  std::endl;
			std::cout << "################################################################  Groesse  " << sweeper->get_end_state()->data().size() <<  std::endl;


#if DIMENSION!=1
    			if(output){
        			GridType::LevelGridView gridView = grid->levelGridView(0);
        			Dune::VTKWriter<GridView> vtkWriter(gridView);
        			Dune::VTKWriter<GridView> vtkWriter2(gridView);
        			string name2 = std::to_string(time);
        			vtkWriter2.addVertexData(sweeper->get_end_state()->data(), "fe_solution_u");
        			vtkWriter2.write("fe_2d_nach_solve" + name2);
			}
#endif

			if((*_new_newton_state[num_nodes]).infinity_norm() < newton){ 

				for(int j=0; j<num_nodes +1 ; j++){
    						_new_initial_state = sweeper->new_newton_state()[j]->data();
    				}
	
				std::cout << "************************************* STARTING NEW TIMESTEP "<< time << std::endl;
	
				break;}



			for(int j=0; j<num_nodes +1 ; j++){
    					(*_new_newton_state[j]) = sweeper->new_newton_state()[j]->data();
    			}
		

    	

   		}//Newtonschritte

	}//Zeitschritte
    
      	auto ut = MPI_Wtime()-st;
        double time;
        MPI_Allreduce(&ut, &time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        std::cout << "Zeit ist am ende " << time << std::endl;

 
}


