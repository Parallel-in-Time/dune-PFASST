#include <config.h>
#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include "../../datatypes/dune_vec.hpp"
#include "fischer_sweeper.hpp"

//#include "assemble.hpp"

#include <dune/fufem/functionspacebases/p1nodalbasis.hh>
#include <dune/fufem/assemblers/operatorassembler.hh>
#include <dune/fufem/assemblers/functionalassembler.hh>
//#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
//#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
//#include <dune/fufem/assemblers/localassemblers/l2functionalassembler.hh>
//#include <dune/fufem/functiontools/basisinterpolator.hh>


using namespace pfasst::examples::fischer_example;

using std::shared_ptr;

using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>; 


using FE_function = Dune::Functions::PQkNodalBasis<GridType::LevelGridView, BASE_ORDER>;  
using sweeper_t = fischer_sweeper<dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>,   FE_function >;

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

    const size_t nelements = get_value<size_t>("num_elements", 16);    // spacial dimension: number of grid points per dimension on the coase level
    
    const double t_0 = 0.0;                                             // left point of the time intervall is zero 
    const double dt = get_value<double>("dt", 0.1);                    // size of timesteping
    double t_end = get_value<double>("tend", 0.1);                      // right point of the time intervall  
    const size_t nnodes = get_value<size_t>("num_nodes", 3);            // time intervall: number of sdc quadrature points
    const QuadratureType quad_type = QuadratureType::GaussRadau;        // quadrature type
    const size_t niter = get_value<size_t>("num_iters", 10);            // maximal number of sdc iterations

    
    
    typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
    typedef GridType::LevelGridView GridView;
    using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
    
    std::shared_ptr<TransferOperatorAssembler<GridType>> transfer;

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
    
    /*auto sdc = std::make_shared<heat_FE_sdc_t>();
	

    
    MatrixType mass;
    MatrixType stiffness;

    auto sweeper = std::make_shared<sweeper_t>(fe_basis[0] , 0, grid); // mass and stiff are just dummies
    sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
    
    sdc->add_sweeper(sweeper);
    sdc->set_options();
    sdc->status()->time() = t_0;
    sdc->status()->dt() = dt;
    sdc->status()->t_end() = t_end;
    sdc->status()->max_iterations() = niter;
    sdc->setup();

    sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());

    //for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << sweeper->initial_state()->data()[i] << std::endl;

    sdc->run();
    
    
    //do not need a post run for GaussRadau nodes not sure about that should ask robert
    //sdc->post_run();

    auto naeherung = sweeper->get_end_state()->data();
    auto exact     = sweeper->exact(t_end)->data();
    auto initial   = sweeper->exact(0)->data();
    for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << initial[i] << " " << naeherung[i] << " " << exact[i] << std::endl;

    sweeper->get_end_state()->scaled_add(-1.0 , sweeper->exact(t_end));
    std::cout << "error in infinity norm: " << sweeper->get_end_state()->norm0()<<  std::endl ;*/


	const auto num_nodes = nnodes;	
    	const auto num_time_steps = 1; //t_end/dt;
	
	vector<vector<shared_ptr<Dune::BlockVector<Dune::FieldVector<double, 1>>>>>  _new_newton_state;
    	_new_newton_state.resize(num_time_steps);

    	for(int i=0; i< num_time_steps; i++){	
		_new_newton_state[i].resize(num_nodes + 1);
		for(int j=0; j<num_nodes +1 ; j++){
			_new_newton_state[i][j] = std::make_shared<Dune::BlockVector<Dune::FieldVector<double, 1>>>(fe_basis[0]->size());
		}
    	}


for(int time=0; time<(t_end-t_0)/dt; time++){	
//int time=0;
    for(int ne=0; ne<4; ne++){	


	auto sweeper = std::make_shared<sweeper_t>(fe_basis[0] , 0, grid); // mass and stiff are just dummies
	sweeper->is_coarse = false;
	std::cout << "test1 "<< std::endl;
    	sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);
	std::cout << "test2 "<< std::endl;
    	auto sdc = std::make_shared<heat_FE_sdc_t>();    
    	sdc->add_sweeper(sweeper);
    	sdc->set_options();
    	sdc->status()->time() = t_0 + time*dt;
    	sdc->status()->dt() = dt;
    	sdc->status()->t_end() = t_0 + (time+1)*dt;
	std::cout << t_0 << " "<< t_0 + time*dt <<" "<< t_0 + (time+1)*dt<< " " << t_end<< std::endl;
    	sdc->status()->max_iterations() = niter;
	std::cout << "vor setup "<< std::endl;
    	sdc->setup();

	std::cout << "test3 "<< std::endl;
	if(time==0 && ne==0) 	//im ersten Newton Lauf Anfangswerte setzen
	for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1; j++){
		for(int k=0; k< _new_newton_state[i][j]->size(); k++){
    		 (*_new_newton_state[i][j])[k]= 0; //sweeper->exact(i*dt+j*1/3.*dt)->data()[k];
    		}
		}
	}
	std::cout << "test4 "<< std::endl;

        sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());

	for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1; j++){
			for(int k=0; k< _new_newton_state[i][j]->size(); k++){
    			sweeper->last_newton_state()[i][j]->data()[k] = (*_new_newton_state[i][j])[k]  ;
			}
    		}
	}
	std::cout << "test5 "<< std::endl;

	    for(int m=0; m< num_nodes +1; m++){
	    	sweeper->df_dune[0][m] = std::make_shared<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>>(sweeper->M_dune); 
            	sweeper->evaluate_df2(*sweeper->df_dune[0][m], sweeper->last_newton_state()[0][m]);
		auto result = sweeper->get_encap_factory().create();
            	result->zero();
                sweeper->evaluate_f2(result, sweeper->last_newton_state()[0][m]);
		sweeper->df_dune[0][m]->mmv(sweeper->last_newton_state()[0][m]->data(), result->data());

	    	sweeper->coarse_rhs()[0][m]->data() =result->data();

	    }


    	sdc->run();   
	sdc->post_run();

	//for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << "+++++++++++++++ new start value " <<sweeper->last_newton_state()[num_time_steps-1][num_nodes ]->data()[i] << " " << (*_new_newton_state[num_time_steps-1][num_nodes])[i]<< " " << sweeper->get_end_state()->data()[i]<< " " << sweeper->states()[num_nodes]->get_data()[i] <<  std::endl;


    	auto naeherung = sweeper->get_end_state()->data();
    	auto exact     = sweeper->exact(sdc->status()->t_end())->data();
    	auto initial   = sweeper->exact(t_0 + time*dt)->data();
    	for(int i=0; i< sweeper->get_end_state()->data().size(); i++) std::cout << initial[i] << " result " << naeherung[i] << " " << exact[i] << std::endl;
	sweeper->get_end_state()->scaled_add(-1.0 , sweeper->exact(sdc->status()->t_end()));
        std::cout << "***************************************    error in infinity norm: " << time << " "<< sweeper->get_end_state()->norm0()<<  std::endl ;

    	for(int i=0; i< num_time_steps; i++){	
		for(int j=0; j<num_nodes +1 ; j++){
		for(int k=0; k< _new_newton_state[i][j]->size(); k++)
    		(*_new_newton_state[i][j])[k] = sweeper->new_newton_state()[i][j]->data()[k];
    		}
	}	

	//std::cout << "################################################################################# this i want to copy " << sweeper->last_newton_state()[num_time_steps-1][num_nodes]->data()[5] << std::endl;

    	

   }
//std::exit(0);
}
    

 
}

