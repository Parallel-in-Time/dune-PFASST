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
#include "../../datatypes/dune_vec.hpp"
#include "../../finite_element_stuff/spectral_transfer.hpp"
//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

//const size_t DIM = 1;            //R??umliche Dimension des Rechengebiets ruth_dim

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

typedef Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, 1, DIMENSION>> SweeperType1;
typedef Heat_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, 2, DIMENSION>> SweeperType2;


typedef pfasst::transfer_traits<SweeperType1, SweeperType2, 1>       TransferTraits;
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

        int nsteps = (t_end-t_0)/dt;
        int max_Newton=2;
          
        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        
        
        
        

        
        auto FinEl = make_shared<fe_manager>(nelements, 2);
        
        vector<vector<shared_ptr<VectorType>>> coarse_all_time_states;
        vector<vector<shared_ptr<VectorType>>> fine_all_time_states;
        
        for(int j=0; j< nsteps; j++){
        coarse_all_time_states[j].resize(nnodes);
        fine_all_time_states[j].resize(nnodes);
        for(int k =0; k< nnodes; j++){
                Dune::BlockVector<Dune::FieldVector<double,NR_COMP> > c =  Dune::BlockVector<Dune::FieldVector<double,NR_COMP> >(FinEl->get_basis1()->size()) ;
                coarse_all_time_states[j][k]=c;
                Dune::BlockVector<Dune::FieldVector<double,NR_COMP> > f =  Dune::BlockVector<Dune::FieldVector<double,NR_COMP> >(FinEl->get_basis2()->size()) ;
                fine_all_time_states[j][k]=f;
        } 
        }
        
        

        for(int i=0; i<max_Newton; i++){
            
            std::cout << "in schleife " << std::endl;
            TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
            pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);
            auto coarse = std::make_shared<SweeperType1>(FinEl->get_basis1(), 1);
            coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
            auto fine = std::make_shared<SweeperType2>(FinEl->get_basis2(), 0);
            fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);
            coarse->is_coarse=true;
            fine->is_coarse=false;
            auto transfer = std::make_shared<TransferType>();
            transfer->create(FinEl);
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
            
            coarse->_old_time_states.resize(nsteps);
            fine->_old_time_states.resize(nsteps);
            for(int j=0; j< nsteps; j++){
                coarse->_old_time_states[j].resize(nnodes);
                fine->_old_time_states[j].resize(nnodes);
                auto& factory = coarse->get_encap_factory();
                std::generate(coarse->_old_time_states[j].begin(), coarse->_old_time_states[j].end(), [&factory](){ return factory.create(); });
                auto& factory2 = fine->get_encap_factory();
                std::generate(fine->_old_time_states[j].begin(), fine->_old_time_states[j].end(), [&factory2](){ return factory2.create(); });
            }
            
            
            
            
            
            for(int j =0; j< nsteps; j++){
               for(int k =0; k< nnodes; j++){
                coarse->_old_time_states[j][k]->data()=_coarse->_old_time_states[j][k]->data();
                fine->_old_time_states[j][k]->data()= _fine->_old_time_states[j][k]->data();
               } 
            }

                        
            
            std::cout << "starte newton iteration " << std::endl;
            
            
            std::cout << "a " << coarse->_old_time_states[0][0]->data()[0] << std::endl;

            std::cout << "vor run " << std::endl;
    
            pfasst.run();
            pfasst.post_run();
            
            std::cout << "a " << coarse->_all_time_states[0][0]->data().size() << std::endl;
            std::cout << "o " << coarse->_old_time_states[0][0]->data().size() << std::endl;


            
            for(int j =0; j< nsteps; j++){
               for(int k =0; k< nnodes; j++){
                _coarse->_old_time_states[j][k]->data()=coarse->_all_time_states[j][k]->data();
                _fine->_old_time_states[j][k]->data()=fine->_all_time_states[j][k]->data();
               } 
            }
            
            //coarse->_all_time_states.clear();
            //fine->_all_time_states.clear();
            
            std::cout << "neues element " << coarse->_old_time_states[0][0]->data()[0] << std::endl;
            
            //coarse->last_newton_state() = coarse->states();
            //fine->last_newton_state() = fine->states();
        
            MPI_Barrier(MPI_COMM_WORLD);
            if(my_rank==num_pro-1) {
                auto naeherung = fine->get_end_state()->data();
                auto exact = fine->exact(t_end)->data();
                for (int i=0; i< fine->get_end_state()->data().size(); i++){
                    std::cout << fine->exact(0)->data()[i] << " " << fine->states()[fine->get_states().size() - 1]->data()[i] << "   " << exact[i] << std::endl;
                }
            std::cout << "******************************************* " << std::endl;
            std::cout << " " << std::endl;
            std::cout << " " << std::endl;
            std::cout << "Fehler: " << std::endl;
            fine->states()[fine->get_states().size() - 1]->scaled_add(-1.0, fine->exact(t_end));
            std::cout << fine->states()[fine->get_states().size() - 1]->norm0() << std::endl;
            ofstream f;
            stringstream ss;
            ss << nelements;
            string s = "solution_pfasst/" + ss.str() + ".dat";
            f.open(s, ios::app | std::ios::out );
            f << nelements << " " << dt << " "<< fine->states()[fine->get_states().size()-1]->norm0() << endl;
            f.close();
            std::cout << "******************************************* " << std::endl;
        
        

        }
            
        }
        
        /*if(my_rank==0) {
          auto naeherung = fine->get_end_state()->data();
          auto exact = fine->exact(t_end)->data();

          /*std::cout << "anzahl states " << fine->get_states().size() << std::endl;
          
          for(int s=0; s< fine->get_states().size(); s++){
              for (int i=0; i< fine->get_end_state()->data().size(); i++){
                std::cout <<  fine->states()[s]->data()[i] << std::endl;
              }
              std::cout <<  "+++++++++++++++++ naester +++++++++++++++++++++" << std::endl;
          }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);*/
        




      }
    }  // ::pfasst::examples::heat_FE
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  MPI_Init(&argc, &argv);

  pfasst::init(argc, argv, SweeperType2::init_opts);
  pfasst::Status<double>::create_mpi_datatype();


  const size_t nelements = get_value<size_t>("num_elements", 4); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("num_nodes", 4);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.1);
  double t_end = get_value<double>("tend", 0.2);
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

  pfasst::examples::heat_FE::run_pfasst(nelements, BASE_ORDER, DIMENSION, nnodes, quad_type, t_0, dt, t_end, niter);

  pfasst::Status<double>::free_mpi_datatype();

  MPI_Finalize();

  return 0;
}
