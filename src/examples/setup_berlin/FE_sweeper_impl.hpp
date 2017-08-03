#include "FE_sweeper.hpp"

#include "assemble.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <utility>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>


using namespace std;
double pi = 3.14159265359;
double two_pi = 2*pi;
double pi_sqr= pi*pi;

namespace pfasst
{
  namespace examples
  {
    namespace FE_sweeper
    {
        
      template<class SweeperTrait, class BasisFunction, typename Enabled>
      void Heat_FE<SweeperTrait, BasisFunction, Enabled>::init_opts(){}
      
      
      template<class SweeperTrait, class BasisFunction, typename Enabled>
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::Heat_FE(std::shared_ptr<BasisFunction> basis, size_t nlevel, std::shared_ptr<GridType> grid) : IMEX<SweeperTrait, BasisFunction, Enabled>()

      {
        std::cout << "fe"<< std::endl;
  
      
	//this->nlevel = nlevel;    
	    
	this->basis = basis;

        //this->grid = grid;
	
        
        //assembleProblem(basis, (this->A_dune), (this->M_dune));
        
        //for(int i =0; i< this->A_dune.N(); i++)
        //for(int j =0; j< this->A_dune.N(); j++)            
        //if (this->A_dune.exists(i,j)) std::cout << this->M_dune[i][j] << " "<< mass[i][j]<<  std::endl;
        //this->A_dune = std::make_shared<MatrixType>(stiffness);
        //this->M_dune = std::make_shared<MatrixType>(mass);
        
        /*std::cout << "vor schleife" <<  std::endl; 
        for(int i =0; i< this->A_dune->N(); i++)
        for(int j =0; j< this->A_dune->N(); j++)            
        //if (this->M_dune->exists(i,j)) std::cout << (*this->M_dune)[i][j] << " "<< mass[i][j] <<std::endl;        
        std::cout << "nach schleife" <<  std::endl; 
        */
        //std::exit(0);
        
        const auto bs = basis->size();
        std::cout << "Finite Element basis of level " << nlevel << " consists of " <<  basis->size() << " elements " << std::endl;

        this->encap_factory()->set_size(bs);

      }

      template<class SweeperTrait, class BasisFunction, typename Enabled>
      void
      Heat_FE<SweeperTrait, BasisFunction,  Enabled>::set_options()
      {
        IMEX<SweeperTrait, BasisFunction, Enabled>::set_options();
        //         this->_nu = config::get_value<spatial_t>("nu", this->_nu);
        //         int num_nodes = this->get_quadrature()->get_num_nodes();
      }

      template<class SweeperTrait, class BasisFunction, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, BasisFunction,  Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();
        std::cout << "need to implement exact " << std::endl;

        
        return result;
      }
      
      template<class SweeperTrait, class BasisFunction, typename Enabled>
      void
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::post_step()
      {
        IMEX<SweeperTrait,BasisFunction,  Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        //this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, class BasisFunction, typename Enabled>
      bool
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, BasisFunction, Enabled>::converged(pre_check);

        if (!pre_check) {
          assert(this->get_status() != nullptr);
          const typename traits::time_t t = this->get_status()->get_time();
          const typename traits::time_t dt = this->get_status()->get_dt();

          assert(this->get_quadrature() != nullptr);
          auto nodes = this->get_quadrature()->get_nodes();
          const auto num_nodes = this->get_quadrature()->get_num_nodes();
          nodes.insert(nodes.begin(), typename traits::time_t(0.0));

          ML_CVLOG(1, this->get_logger_id(),
                   "Observables after "
                   << ((this->get_status()->get_iteration() == 0)
                          ? std::string("prediction")
                          : std::string("iteration ") + std::to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(),
                     "  t["<<m<<"]=" <<  (t + dt * nodes[m])
                     << "      |abs residual| = " <<  this->_abs_res_norms[m]
                     << "      |rel residual| = " <<  this->_rel_res_norms[m]
                    //                      << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
                    //                      << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" <<  (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " <<  this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " <<  this->_rel_res_norms[num_nodes]
//                   << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
//                   << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes])
                 );
        }
        return converged;
      }

      template<class SweeperTrait, class BasisFunction, typename Enabled>
      bool
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, class BasisFunction, typename Enabled>
      size_t
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }

      template<class SweeperTrait, class BasisFunction, typename Enabled> // calculates the error on every quadrature point
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
      {
        ML_CVLOG(4, this->get_logger_id(), "computing error");

        assert(this->get_status() != nullptr);
        const typename traits::time_t dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), typename traits::time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> error;
        error.resize(num_nodes + 1);
        std::generate(error.begin(), error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const typename traits::time_t ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), this->get_states()[m]);
        }

        return error;
      }

      template<class SweeperTrait, class BasisFunction, typename Enabled> 
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                            const typename SweeperTrait::time_t& t)
      {
        UNUSED(t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> rel_error;
        rel_error.resize(error.size());
        std::generate(rel_error.begin(), rel_error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          rel_error[m]->scaled_add(1.0 / this->get_states()[m]->norm0(), error[m]);
        }

        return rel_error;
      }

      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t> Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t, const shared_ptr<typename SweeperTrait::encap_t> u)      { }*/

      template<class SweeperTrait, class BasisFunction, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
	
        ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);
        auto result = this->get_encap_factory().create();        
        return result;
      }

      
      
      
      template<class SweeperTrait, class BasisFunction, typename Enabled>
      void
      Heat_FE<SweeperTrait, BasisFunction, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {ML_CVLOG(4, this->get_logger_id(), "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);  }
      

      
      
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
