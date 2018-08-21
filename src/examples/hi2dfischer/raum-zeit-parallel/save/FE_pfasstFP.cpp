//#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
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


#include "FE_sweeper.hpp"
#include <pfasst/encap/dune_vec.hpp>
#include "../1d_transfer/spectral_transfer.hpp"
//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

//const size_t DIM = 1;            //Räumliche Dimension des Rechengebiets ruth_dim

//const size_t BASIS_ORDER = 1;    //maximale Ordnung der Lagrange Basisfunktionen

//////////////////////////////////////////////////////////////////////////////////////


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;
using pfasst::encap::DuneEncapsulation;
//using pfasst::encap::VectorEncapsulation;

using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::contrib::SpectralTransfer;
using pfasst::TwoLevelPfasst;
typedef pfasst::comm::MpiP2P CommunicatorType;

using pfasst::examples::heat_FE::Heat_FE;

typedef DuneEncapsulation<double, double, 1>                     EncapType;
//typedef VectorEncapsulation<double, double, 1>                     EncapType;


//typedef Heat_FE<pfasst::sweeper_traits<typename EncapType::traits>> SweeperType;
typedef Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, BASE_ORDER, DIMENSION>> SweeperType;


typedef pfasst::transfer_traits<SweeperType, SweeperType, 1>       TransferTraits;
typedef SpectralTransfer<TransferTraits>                           TransferType;


namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      void run_pfasst(const size_t nelements, const size_t basisorder, const size_t dim, const size_t& nnodes, const pfasst::quadrature::QuadratureType& quad_type,
                      const double& t_0, const double& dt, const double& t_end, const size_t& niter, double newton)
      {


        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        
        MPI_Comm comm_x;//=MPI_COMM_SELF, 
        MPI_Comm comm_t;//=MPI_COMM_WORLD; 
	int myid, xcolor, tcolor;

	int space_num=2;
	int space_size, time_size;
   	xcolor = (my_rank / space_num);
   	tcolor = my_rank % space_num;



        
        //rank   xcolor    tcolor
        // 0       0         0
        // 1       0         1
        // 2       1         0
        // 3       1         1

   	MPI_Comm_split( MPI_COMM_WORLD, xcolor, my_rank, &comm_x );
   	MPI_Comm_split( MPI_COMM_WORLD, tcolor, my_rank, &comm_t );
	
	int space_rank=77;
        MPI_Comm_rank(comm_x, &space_rank );
        MPI_Comm_size(comm_x, &space_size);
        
        int time_rank=77;
        MPI_Comm_rank(comm_t, &time_rank );	

	
	std::cout << my_rank << space_rank << time_rank << " space_size "<<space_size << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	//std::exit(0);
	


        TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
        pfasst.communicator() = std::make_shared<CommunicatorType>(comm_t);
        
        //pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);
        //pfasst.grid_builder(nelements);
	auto FinEl = make_shared<fe_manager>(nelements, 2, comm_x);

        auto coarse = std::make_shared<SweeperType>(FinEl, 1);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<SweeperType>(FinEl, 0);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        auto transfer = std::make_shared<TransferType>();
	transfer->create(FinEl);

        coarse->is_coarse=true;
        fine->is_coarse=false;
        coarse->comm=comm_x;
        fine->comm=comm_x;
    
    
        //pfasst.add_sweeper(coarse, true);
        //pfasst.add_sweeper(fine, false);
        pfasst.add_transfer(transfer);
	pfasst.add_sweeper(coarse, true);
	pfasst.add_sweeper(fine);
        pfasst.set_options();

        //std::cout << "hier im code" <<std::endl;

        pfasst.status()->time() = t_0;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_end;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

        coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());
        fine->initial_state() = fine->exact(pfasst.get_status()->get_time());
	fine->newton=newton;
	coarse->newton=newton;



    	MPI_Barrier(MPI_COMM_WORLD);
    	auto st = MPI_Wtime();
        pfasst.run();
        pfasst.post_run();
    	auto ut = MPI_Wtime()-st;
    	double time;
    	MPI_Reduce(&ut, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    	std::cout << my_rank << "Zeit ist am end" << time << std::endl;



	int global_num;
	MPI_Reduce(&(fine->num_solves), &global_num, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


        //if(space_rank==space_size-1) {

          //std::cout << "******************************************* " << std::endl;
          //std::cout << " " << std::endl;

          auto naeherung  = fine->get_end_state()->data();
          auto exact      = fine->exact(t_end)->data();
          auto initial    = fine->exact(0)->data();          
          /*for (int i=0; i< fine->get_end_state()->data().size(); i++){
            if(time_rank==1 && space_rank==1) std::cout << initial[i] << " " << naeherung[i] << "   " << exact[i] << std::endl;
          }*/
          /*for (int i=0; i< sweeper->get_end_state()->data().size(); i++){
            if(rank==1) std::cout << sweeper->exact(0)->data()[i] << " " << naeherung[i] << "   " << exact[i] << std::endl;
          }*/
          //MPI_Barrier(comm_t);
          //std::cout << " " << std::endl;
          //std::cout << "Fehler: " << std::endl;
          //auto norm =  fine->exact(t_end))->data();
          fine->end_state()->scaled_add(-1.0, fine->exact(t_end));
          /*MPI_Barrier(MPI_COMM_WORLD);          
          for (int i=0; i< fine->get_end_state()->data().size(); i++){
            if(time_rank==1 && space_rank==1) std::cout << fine->end_state()->get_data()[i] << std::endl;
          }
          MPI_Barrier(MPI_COMM_WORLD);*/
          //std::cout << time_rank << " " << space_rank << " der lokale fehler betraegt " << fine->end_state()->norm0() << std::endl;          
          
          
          double error= fine->end_state()->norm0(true, comm_x);
          std::cout << time_rank << " " << space_rank << " der fehler betraegt " << error << std::endl;
        
        //}

	std::cout << my_rank << " num_solves " << fine->num_solves << std::endl;


        MPI_Barrier(MPI_COMM_WORLD);

        /*for (int i=0; i<num_pro; i++){
          if(my_rank==i){
            ofstream ff;
            stringstream sss;
            sss << nelements << "_iter";
            string st = "solution_pfasst/" + sss.str() + ".dat";
            ff.open(st, ios::app | std::ios::out );
            auto iter = pfasst._it_per_step;
            for (const auto &line : iter) {
              ff << my_rank << " " << dt <<"     " << line << std::endl;

            }

            ff.close();

          }
          MPI_Barrier(MPI_COMM_WORLD);

        }*/






      }
    }  // ::pfasst::examples::heat_FE
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  MPI_Init(&argc, &argv);

  pfasst::init(argc, argv, SweeperType::init_opts);
  pfasst::Status<double>::create_mpi_datatype();


  const size_t nelements = get_value<size_t>("num_elements", 4); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0;
  const double dt = get_value<double>("dt", 0.1);
  double t_end = get_value<double>("tend", 0.2);
  size_t nsteps = get_value<size_t>("num_steps", 0);
  double newton = get_value<double>("newton", 1e-2);
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
  const size_t niter = get_value<size_t>("num_iters", 20);







  pfasst::examples::heat_FE::run_pfasst(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter, newton);



  pfasst::Status<double>::free_mpi_datatype();

  MPI_Finalize();

  return 0;
}
