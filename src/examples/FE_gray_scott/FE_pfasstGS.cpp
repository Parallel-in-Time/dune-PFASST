//#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

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
#include "dune_vec_multi.hpp"
#include "../../finite_element_stuff/spectral_transfer_gs.hpp"


//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

const size_t DIM = 2;            //Räumliche Dimension des Rechengebiets ruth_dim

const size_t BASIS_ORDER = 1;    //maximale Ordnung der Lagrange Basisfunktionen

//////////////////////////////////////////////////////////////////////////////////////


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 2, 2>;  //

//using pfasst::encap::DuneEncapsulation;


using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::contrib::SpectralTransfer;
using pfasst::TwoLevelPfasst;
typedef pfasst::comm::MpiP2P CommunicatorType;

using pfasst::examples::heat_FE::Heat_FE;

//typedef DuneEncapsulation<double, double, 2,2>  EncapType;




typedef Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM, 2>> SweeperType;


typedef pfasst::transfer_traits<SweeperType, SweeperType, 1>       TransferTraits;
typedef SpectralTransfer<TransferTraits>                           TransferType;


namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      void run_pfasst(const size_t nelements, const size_t basisorder, const size_t dim, const size_t& nnodes, const pfasst::quadrature::QuadratureType& quad_type,
                      const double& t_0, const double& dt, const double& t_end, const size_t& niter)
      {


        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
        pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);
        

	auto FinEl = make_shared<fe_manager>(nelements, 2);





        auto coarse = std::make_shared<SweeperType>(FinEl, 1);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        std::cout << "grober sweeper erstellt : " << my_rank<< std::endl;

        auto fine = std::make_shared<SweeperType>(FinEl, 0);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        std::cout << "feiner sweeper erstellt : " << my_rank<< std::endl;

        auto transfer = std::make_shared<TransferType>();

        transfer->create(FinEl);

        std::cout << "2 Hallo Welt ich bin Prozess : " << my_rank<< std::endl;

        pfasst.add_sweeper(coarse, true);
        pfasst.add_sweeper(fine);


        pfasst.add_transfer(transfer);
        std::cout << "3 Hallo Welt ich bin Prozess : " << my_rank<< std::endl;
        //pfasst.add_sweeper(coarse, fine);
        std::cout << "4 Hallo Welt ich bin Prozess : " << my_rank<< std::endl;
        pfasst.set_options();



        std::cout << "5 Hallo Welt ich bin Prozess : " << my_rank<< std::endl;

        pfasst.status()->time() = t_0;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_end;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

        coarse->is_coarse= true;
        fine->is_coarse=false;

        coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());
        fine->initial_state() = fine->exact(pfasst.get_status()->get_time());

        //std::cout << "vor run Hallo Welt ich bin Prozess : " << my_rank<< std::endl;

	double t1, t2;
	t1=MPI_Wtime();	
        
	pfasst.run();
        pfasst.post_run();

	t2=MPI_Wtime();
	printf("time is %f\n", t2-t1);


        MPI_Barrier(MPI_COMM_WORLD);

        for (int i=0; i<num_pro; i++){
          if(my_rank==i){
            ofstream ff;
            stringstream sss;
	    sss << num_pro << "_" << nelements*2 << "_" << dt << "_" << t_end;	
            string st = "solution_pfasst/" + sss.str() + ".dat";
            ff.open(st, ios::app | std::ios::out );
            auto iter = pfasst._it_per_step;
            for (const auto &line : iter) {
              ff << my_rank << " " << line << std::endl;
            }
 	    ff << "  time = " << time1 << endl;
            ff.close();

          }
          MPI_Barrier(MPI_COMM_WORLD);

        }



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


  const size_t nelements = get_value<size_t>("num_elements", 32); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 2);
  double t_end = get_value<double>("tend", 12);
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
  const size_t niter = get_value<size_t>("num_iters", 50);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
  std::cout << "Hallo Welt ich bin Prozess : " << my_rank<< std::endl;

  pfasst::examples::heat_FE::run_pfasst(nelements, BASIS_ORDER, DIM, nnodes, quad_type, t_0, dt, t_end, niter);

  pfasst::Status<double>::free_mpi_datatype();

  MPI_Finalize();

  return 0;
}
