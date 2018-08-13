
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


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;
using pfasst::encap::DuneEncapsulation;


using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::contrib::SpectralTransfer;
using pfasst::TwoLevelPfasst;
typedef pfasst::comm::MpiP2P CommunicatorType;

using pfasst::examples::heat_FE::Heat_FE;

typedef DuneEncapsulation<double, double, 1>                     EncapType;



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
                      const double& t_0, const double& dt, const double& t_end, const size_t& niter, double newton, bool output)
      {


        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
        pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);
	auto FinEl = make_shared<fe_manager>(nelements, 2);

        auto coarse = std::make_shared<SweeperType>(FinEl, 1);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<SweeperType>(FinEl, 0);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        auto transfer = std::make_shared<TransferType>();
	transfer->create(FinEl);

        coarse->is_coarse=true;
        fine->is_coarse=false;

        pfasst.add_transfer(transfer);
	pfasst.add_sweeper(coarse, true);
	pfasst.add_sweeper(fine);
        pfasst.set_options();

        pfasst.status()->time() = t_0;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_end;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

        coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());
        fine->initial_state() = fine->exact(pfasst.get_status()->get_time());
	fine->newton=newton;
	coarse->newton=newton;

        pfasst.run();
        pfasst.post_run();

	int global_num;
	MPI_Reduce(&(fine->num_solves), &global_num, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


        if(output && my_rank==num_pro-1) {
          auto naeherung = fine->get_end_state()->data();
          auto exact = fine->exact(t_end)->data();
          for (int i = 0; i < fine->get_end_state()->data().size(); i++) {
            std::cout << fine->exact(0)->data()[i] << " " << naeherung[i] << "   " << exact[i] << std::endl;
          }
          std::cout << "******************************************* " << std::endl;
          std::cout << " " << std::endl;
          std::cout << " " << std::endl;
        }
        std::cout << "ERROR (infinity norm): " << std::endl;
        fine->states()[fine->get_states().size() - 1]->scaled_add(-1.0, fine->exact(t_end));
        std::cout << fine->states()[fine->get_states().size() - 1]->norm0() << std::endl;
	std::cout << "my_rank: " << my_rank << ", number solutions which this process does: " << fine->num_solves << "" << std::endl;


        MPI_Barrier(MPI_COMM_WORLD);

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
  bool output = get_value<double>("output", 0);
  const size_t niter = get_value<size_t>("num_iters", 1000);

  pfasst::examples::heat_FE::run_pfasst(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter, newton, output);

  pfasst::Status<double>::free_mpi_datatype();

  MPI_Finalize();

  return 0;
}
