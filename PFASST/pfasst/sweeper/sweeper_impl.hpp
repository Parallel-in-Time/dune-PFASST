#include "pfasst/sweeper/sweeper.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"
using pfasst::quadrature::IQuadrature;


namespace pfasst
{
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  Sweeper<SweeperTrait, BaseFunction, Enabled>::Sweeper()
    :   _quadrature(nullptr)
      , _factory(std::make_shared<typename SweeperTrait::encap_t::factory_t>())
      , _states(0)
      , _previous_states(0)
      , _M_initial(nullptr)
      , _end_state(nullptr)
      , _tau(0)
      , _residuals(0)
      , _status(nullptr)
      , _abs_residual_tol(0.0)
      , _rel_residual_tol(0.0)
      , _logger_id("SWEEPER")
  {}

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<IQuadrature<typename SweeperTrait::time_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::quadrature()
  {
    return this->_quadrature;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const shared_ptr<IQuadrature<typename SweeperTrait::time_t>>
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_quadrature() const
  {
    return this->_quadrature;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<Status<typename SweeperTrait::time_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::status()
  {
    return this->_status;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const shared_ptr<Status<typename SweeperTrait::time_t>>
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_status() const
  {
    return this->_status;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_t::factory_t>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::encap_factory()
  {
    return this->_factory;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const typename SweeperTrait::encap_t::factory_t&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_encap_factory() const
  {
    return *(this->_factory);
  }

  /**
   * @throws std::runtime_error if `get_states()` has zero length
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_t>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::initial_state()
  {
    if (this->get_states().size() == 0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Sweeper need to be setup first before querying initial state.");
      throw std::runtime_error("sweeper not setup before querying initial state");
    }
    return this->states().front();
  }

  /**
   * @throws std::runtime_error if `get_states()` has zero length
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const shared_ptr<typename SweeperTrait::encap_t>
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_initial_state() const
  {
    if (this->get_states().size() == 0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Sweeper need to be setup first before querying initial state.");
      throw std::runtime_error("sweeper not setup before querying initial state");
    }
    return this->get_states().front();
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::states()
  {
    return this->_states;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_states() const
  {
    return this->_states;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::previous_states()
  {
    return this->_previous_states;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_previous_states() const
  {
    return this->_previous_states;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_t>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::end_state()
  {
    return this->_end_state;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const shared_ptr<typename SweeperTrait::encap_t>
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_end_state() const
  {
    return this->_end_state;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::tau()
  {
    return this->_tau;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_tau() const
  {
    return this->_tau;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::residuals()
  {
    return this->_residuals;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_t>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_residuals() const
  {
    return this->_residuals;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::set_logger_id(const string& logger_id)
  {
    this->_logger_id = logger_id;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const char*
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_logger_id() const
  {
    return this->_logger_id.c_str();
  }

  /**
   * @note Sets tolerances for absolute and relative residual if given on the command line.
   *   Otherwise the currently set values are unchanged.
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::set_options()
  {
    ML_CVLOG(3, this->get_logger_id(), "setting options from runtime parameters (if available)");
    this->_abs_residual_tol = config::get_value<typename traits::spatial_t>("abs_res_tol", this->_abs_residual_tol);
    this->_rel_residual_tol = config::get_value<typename traits::spatial_t>("rel_res_tol", this->_rel_residual_tol);
    ML_CVLOG(3, this->get_logger_id(), "  absolute residual tolerance: " << this->_abs_residual_tol);
    ML_CVLOG(3, this->get_logger_id(), "  relative residual tolerance: " << this->_rel_residual_tol);
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction,Enabled>::set_abs_residual_tol(const typename SweeperTrait::spatial_t& abs_res_tol)
  {
    this->_abs_residual_tol = abs_res_tol;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::set_rel_residual_tol(const typename SweeperTrait::spatial_t& rel_res_tol)
  {
    this->_rel_residual_tol = rel_res_tol;
  }

  /**
   * @throws std::runtime_error if either `get_status()` or `get_quadrature()` are not set, i.e. `nullptr`.
   */
  template<class SweeperTrait, class BaseFunction,  typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::initialize()
  {
    if (this->get_status() == nullptr) {
      //throw std::runtime_error("Status not yet set.");
    }
    ML_CVLOG(1, this->get_logger_id(), "setting up with t0=" << this->get_status()->get_time()
              << ", dt=" << this->get_status()->get_dt()
              << ", t_end=" << this->get_status()->get_t_end()
              << ", max_iter=" << this->get_status()->get_max_iterations());

    if (this->get_quadrature() == nullptr) {
      throw std::runtime_error("Quadrature not yet set.");
    }
    ML_CLOG(INFO, this->get_logger_id(), "using as quadrature: " << this->get_quadrature()->print_summary()
                                      //<< " and an expected error of " << LOG_FLOAT << this->get_quadrature()->expected_error());
                                        << " and an expected error of " << this->get_quadrature()->expected_error());



    const auto nodes = this->get_quadrature()->get_nodes();
    const auto num_nodes = this->get_quadrature()->get_num_nodes();

    
    
    this->states().resize(num_nodes + 1);
    auto& factory = this->get_encap_factory();
    //factory.create();
    std::generate(this->states().begin(), this->states().end(), [&factory](){ return factory.create(); });

    this->previous_states().resize(num_nodes + 1);
    std::generate(this->previous_states().begin(), this->previous_states().end(), [&factory](){ return factory.create(); });


    _M_initial= this->get_encap_factory().create();
    
    this->end_state() = this->get_encap_factory().create();

    this->tau().resize(num_nodes + 1);
    std::generate(this->tau().begin(), this->tau().end(), [&factory](){ return factory.create(); });

    this->residuals().resize(num_nodes + 1);
    std::generate(this->residuals().begin(), this->residuals().end(), [&factory](){ return factory.create(); });
  }

  /**
   * @throws std::runtime_error if either `get_status()` or `get_quadrature()` are not set, i.e. `nullptr`.
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::setup()
  {
    if (this->get_status() == nullptr) {
      //throw std::runtime_error("Status not yet set.");
    }
    ML_CVLOG(1, this->get_logger_id(), "setting up with t0=" << this->get_status()->get_time()
              << ", dt=" << this->get_status()->get_dt()
              << ", t_end=" << this->get_status()->get_t_end()
              << ", max_iter=" << this->get_status()->get_max_iterations());

    if (this->get_quadrature() == nullptr) {
      throw std::runtime_error("Quadrature not yet set.");
    }
    ML_CLOG(INFO, this->get_logger_id(), "using as quadrature: " << this->get_quadrature()->print_summary()
                                      //<< " and an expected error of " << LOG_FLOAT << this->get_quadrature()->expected_error());
                                      << " and an expected error of " << this->get_quadrature()->expected_error());

    //this->initialize();
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::reset()
  {
    ML_CLOG(WARNING, this->get_logger_id(),
            "Resetting to emulate recovering from faulty process.");

    this->status()->reset();
    this->initialize();
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::pre_predict()
  {
    ML_CVLOG(4, this->get_logger_id(), "pre-predicting");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::predict()
  {
    ML_CVLOG(4, this->get_logger_id(), "predicting");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::post_predict()
  {
    ML_CVLOG(4, this->get_logger_id(), "post-predicting");

    assert(this->get_status() != nullptr);
    this->integrate_end_state(this->get_status()->get_dt());

//    assert(this->get_quadrature() != nullptr);
//     ML_CVLOG(2, this->get_logger_id(), "solution at nodes:");
//     for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
//       ML_CVLOG(2, this->get_logger_id(), "  " << m << ": " << to_string(this->get_states()[m]));
//     }
//     ML_CVLOG(1, this->get_logger_id(), "solution at t_end: " << to_string(this->get_end_state()));
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::pre_sweep()
  {
    ML_CVLOG(4, this->get_logger_id(), "pre-sweeping");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::sweep()
  {
    ML_CVLOG(4, this->get_logger_id(), "sweeping");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::post_sweep()
  {
    ML_CVLOG(4, this->get_logger_id(), "post-sweeping");

    assert(this->get_status() != nullptr);
    this->integrate_end_state(this->get_status()->get_dt());

//    assert(this->get_quadrature() != nullptr);
//     ML_CVLOG(2, this->get_logger_id(), "solution at nodes:");
//     for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
//       ML_CVLOG(2, this->get_logger_id(), "\t" << m << ": " << to_string(this->get_states()[m]));
//     }
//     ML_CVLOG(1, this->get_logger_id(), "solution at t_end:" << to_string(this->get_end_state()));
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::post_step()
  {
    ML_CVLOG(4, this->get_logger_id(), "post step");

//    assert(this->get_quadrature() != nullptr);
//     ML_CVLOG(2, this->get_logger_id(), "solution at nodes:");
//     for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
//       ML_CVLOG(2, this->get_logger_id(), "\t" << m << ": " << to_string(this->get_states()[m]));
//     }
//     ML_CVLOG(1, this->get_logger_id(), "solution at t_end: " << to_string(this->get_end_state()));
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::advance(const size_t& num_steps)
  {
    ML_CVLOG(1, this->get_logger_id(), "advancing " << num_steps << " time steps");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::advance()
  {
    this->advance(1);
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::spread()
  {
    ML_CVLOG(4, this->get_logger_id(), "spreading initial value to all states");

    assert(this->get_initial_state() != nullptr);
    //std::cout << "anzahl der states "<< this->get_states().size() << std::endl; std::exit(0);	
    for(size_t m = 1; m < this->get_states().size(); ++m) {
      assert(this->states()[m] != nullptr);
      this->states()[m]->data() = this->get_initial_state()->get_data();
    }
  }


  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::spread(double value)
  {
    ML_CVLOG(4, this->get_logger_id(), "spreading initial value to all states");

    assert(this->get_initial_state() != nullptr);
    std::cout << "vor der schleife "<< this->get_states().size() << std::endl; 	
    for(size_t m = 1; m < this->get_states().size(); ++m) {
      assert(this->states()[m] != nullptr);
      this->states()[m]->data() = value; //this->get_initial_state()->get_data();
    }
    std::cout << "nach der schleife "<< this->get_states().size() << std::endl; 	
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::spread_Newton()
  {
    ML_CVLOG(4, this->get_logger_id(), "spreading initial value to all states");

    assert(this->get_initial_state() != nullptr);
    //std::cout << "anzahl der states "<< this->get_states().size() << std::endl; std::exit(0);	
    std::cout << "vor der schleife "<< this->get_states().size() << std::endl; 	
    int i= 0; //(this->get_status()->get_time()/this->get_status()->get_dt()); 	
    for(size_t m = 1; m < this->get_states().size(); ++m) {
      assert(this->states()[m] != nullptr);
      //std::cout << i << " " << this->last_newton_state()[i][m]->data()  << " " << i << std::endl;
      this->states()[m]->data() = this->last_newton_state()[i][m]->data();
    }
    std::cout << "nach der schleife "<< this->get_states().size() << std::endl; 	
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<vector<shared_ptr<typename SweeperTrait::encap_t>>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::last_newton_state()
  {
    return this->_last_newton_state;
  }
  

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const vector<vector<shared_ptr<typename SweeperTrait::encap_t>>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_last_newton_state() const
  {
    return this->_last_newton_state;
  }


  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<vector<shared_ptr<typename SweeperTrait::encap_t>>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::new_newton_state()
  {
    return this->_new_newton_state;
  }
  

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const vector<vector<shared_ptr<typename SweeperTrait::encap_t>>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_new_newton_state() const
  {
    return this->_new_newton_state;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<vector<shared_ptr<typename SweeperTrait::encap_t>>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::coarse_rhs()
  {
    return this->_coarse_rhs;
  }
  //

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  const vector<vector<shared_ptr<typename SweeperTrait::encap_t>>>&
  Sweeper<SweeperTrait, BaseFunction, Enabled>::get_coarse_rhs() const
  {
    return this->_coarse_rhs;
  }


  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::save()
  {
    ML_CVLOG(4, this->get_logger_id(), "saving states to previous states");

    assert(this->get_quadrature() != nullptr);
    assert(this->get_states().size() == this->get_quadrature()->get_num_nodes() + 1);
    assert(this->get_previous_states().size() == this->get_states().size());

       
    _all_time_states.push_back(this->get_states());
    
    
    for (size_t m = 0; m < this->get_states().size(); ++m) {
      this->previous_states()[m]->data() = this->get_states()[m]->get_data();
    }
  }

  /**
   * @throws std::runtime_error if not overwritten in specialized implementation
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::reevaluate(const bool initial_only)
  {
    UNUSED(initial_only);
    throw std::runtime_error("reevaluation of right-hand-side");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::reevaluate()
  {
    this->reevaluate(false);
  }

  /**
   * @throws std::runtime_error if not overwritten in specialized implementation
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>
  Sweeper<SweeperTrait, BaseFunction, Enabled>::integrate(const typename SweeperTrait::time_t& dt)
  {
    UNUSED(dt);
    throw std::runtime_error("integration over dt");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  bool
  Sweeper<SweeperTrait, BaseFunction, Enabled>::converged(const bool pre_check)
  {

    this->compute_residuals(pre_check);
    const size_t num_residuals = this->get_residuals().size();
    this->_abs_res_norms.resize(num_residuals);
    this->_rel_res_norms.resize(num_residuals);


    assert(this->get_residuals().back() != nullptr);

    this->_abs_res_norms.back() = this->get_residuals().back()->norm0();
    //std::cout <<"hier abs res " << this->get_states().back()->norm0() << std::endl;
    this->_rel_res_norms.back() = this->_abs_res_norms.back() / this->get_states().back()->norm0();
    if (pre_check) {
      if (this->_abs_residual_tol > 0.0 || this->_rel_residual_tol > 0.0) {
        ML_CVLOG(4, this->get_logger_id(), "preliminary convergence check");

        if (this->_abs_res_norms.back() < this->_abs_residual_tol) {
          ML_CVLOG(1, this->get_logger_id(), "Sweeper has converged w.r.t. absolute residual tolerance: 1 "
                                          << this->_abs_res_norms.back() << " < " << this->_abs_residual_tol);
        } else if (this->_rel_res_norms.back() < this->_rel_residual_tol) {
          ML_CVLOG(1, this->get_logger_id(), "Sweeper has converged w.r.t. relative residual tolerance: 2 "
                                          << this->_rel_res_norms.back() << " < " << this->_rel_residual_tol);
        } else {
          ML_CVLOG(1, this->get_logger_id(), "Sweeper has not yet converged to neither residual tolerance. 3 ");
        }

        return (   this->_abs_res_norms.back() < this->_abs_residual_tol
                || this->_rel_res_norms.back() < this->_rel_residual_tol);
      } else {
        ML_CLOG(WARNING, this->get_logger_id(), "No residual tolerances set. Thus skipping convergence check.");
        return false;
      }
    } else {
      for (size_t m = 0; m < num_residuals - 1; ++m) {
        assert(this->get_residuals()[m] != nullptr);
        const auto norm = this->get_residuals()[m]->norm0();
        this->_abs_res_norms[m] = norm;
        this->_rel_res_norms[m] = this->_abs_res_norms[m] / this->get_states()[m]->norm0();
      }

      this->status()->abs_res_norm() = *(std::max_element(this->_abs_res_norms.cbegin(), this->_abs_res_norms.cend()));
      this->status()->rel_res_norm() = *(std::max_element(this->_rel_res_norms.cbegin(), this->_rel_res_norms.cend()));

      if (this->_abs_residual_tol > 0.0 || this->_rel_residual_tol > 0.0) {
        ML_CVLOG(4, this->get_logger_id(), "convergence check");

        if (this->status()->abs_res_norm() < this->_abs_residual_tol) {
          ML_CLOG(INFO, this->get_logger_id(), "Sweeper has converged w.r.t. absolute residual tolerance: 4 "
                                          << this->status()->abs_res_norm() << " < " << this->_abs_residual_tol);
        } else if (this->status()->rel_res_norm() < this->_rel_residual_tol) {
          ML_CLOG(INFO, this->get_logger_id(), "Sweeper has converged w.r.t. relative residual tolerance: 5 "
                                          << this->status()->rel_res_norm() << " < " << this->_rel_residual_tol);
        } else {
          ML_CLOG(INFO, this->get_logger_id(), "Sweeper has not yet converged to neither residual tolerance.");
        }

        return (   this->status()->abs_res_norm() < this->_abs_residual_tol
                || this->status()->rel_res_norm() < this->_rel_residual_tol);
      } else {
        ML_CLOG(WARNING, this->get_logger_id(), "No residual tolerances set. Thus skipping convergence check.");
        return false;
      }
    }
  }


    /*template<class SweeperTrait, class BaseFunction, typename Enabled>
    bool
    Sweeper<SweeperTrait, BaseFunction, Enabled>::alternative_converged(const bool pre_check) //u^{k+1} - u^k, wird mit false aufgerufen
    {
      this->compute_residuals(pre_check);

      const size_t num_residuals = this->get_residuals().size();
      this->_abs_res_norms.resize(num_residuals);
      this->_rel_res_norms.resize(num_residuals);

      assert(this->get_residuals().back() != nullptr);
      this->_abs_res_norms.back() = this->get_residuals().back()->norm0();
      this->_rel_res_norms.back() = this->_abs_res_norms.back() / this->get_states().back()->norm0();

      if (pre_check) {
        if (this->_abs_residual_tol > 0.0 || this->_rel_residual_tol > 0.0) {
          ML_CVLOG(4, this->get_logger_id(), "preliminary convergence check");

          if (this->_abs_res_norms.back() < this->_abs_residual_tol) {
            ML_CVLOG(1, this->get_logger_id(), "Sweeper has converged w.r.t. absolute residual tolerance: "
                                               << this->_abs_res_norms.back() << " < " << this->_abs_residual_tol);
          } else if (this->_rel_res_norms.back() < this->_rel_residual_tol) {
            ML_CVLOG(1, this->get_logger_id(), "Sweeper has converged w.r.t. relative residual tolerance: "
                                               << this->_rel_res_norms.back() << " < " << this->_rel_residual_tol);
          } else {
            ML_CVLOG(1, this->get_logger_id(), "Sweeper has not yet converged to neither residual tolerance.");
          }

          return (   this->_abs_res_norms.back() < this->_abs_residual_tol
                     || this->_rel_res_norms.back() < this->_rel_residual_tol);
        } else {
          ML_CLOG(WARNING, this->get_logger_id(), "No residual tolerances set. Thus skipping convergence check.");
          return false;
        }
      } else {
        for (size_t m = 0; m < num_residuals - 1; ++m) {
          assert(this->get_residuals()[m] != nullptr);
          const auto norm = this->get_residuals()[m]->norm0();
          this->_abs_res_norms[m] = norm;
          this->_rel_res_norms[m] = this->_abs_res_norms[m] / this->get_states()[m]->norm0();
        }

        this->status()->abs_res_norm() = *(std::max_element(this->_abs_res_norms.cbegin(), this->_abs_res_norms.cend()));
        this->status()->rel_res_norm() = *(std::max_element(this->_rel_res_norms.cbegin(), this->_rel_res_norms.cend()));
        /////////////////


        vector<typename traits::spatial_t>  new_abs_res_norm(num_residuals-1);
        vector<typename traits::spatial_t>  new_rel_res_norm(num_residuals-1);



        //std::cout << "  abs res norm " << this->status()->abs_res_norm()<< std::endl;
        for (size_t m = 0; m < num_residuals - 1 ; ++m) {

          //auto new_abs_res = &(this->get_previous_states()[m]);
          typename traits::encap_t new_abs_res;
          new_abs_res.data() = this->get_previous_states()[m]->data();
          //SweeperTrait::encap_t new_abs_res = new SweeperTrait::encap_t(this->get_previous_states()[m]->data());
          //std::cout << this->get_states()[m]->norm0() << " " << this->get_previous_states()[m]->norm0() << std::endl;
          new_abs_res.scaled_add(-1., this->get_states()[m]);

          new_abs_res_norm[m]= new_abs_res.norm0();
          new_rel_res_norm[m]= 0; //new_abs_res_norm[m]    / this->get_states()[m]->norm0();
        }
        //std::cout << this->get_previous_states()[0]->norm0() << " norm " << std::endl;

        auto new_abs_res_norm_max = *(std::max_element(new_abs_res_norm.cbegin(), new_abs_res_norm.cend()));
        auto new_rel_res_norm_max = *(std::max_element(new_rel_res_norm.cbegin(), new_rel_res_norm.cend()));

        //this->status()->abs_res_norm() = *(std::max_element(this->_abs_res_norms.cbegin(), this->_abs_res_norms.cend()));
        //this->status()->rel_res_norm() = *(std::max_element(this->_rel_res_norms.cbegin(), this->_rel_res_norms.cend()));

        if (this->_abs_residual_tol > 0.0 || this->_rel_residual_tol > 0.0) {
          ML_CVLOG(4, this->get_logger_id(), "convergence check");

          if (new_abs_res_norm_max < this->_abs_residual_tol) {
            ML_CLOG(INFO, this->get_logger_id(), " uk+1- uk Sweeper has converged w.r.t. absolute residual tolerance: "   << new_abs_res_norm_max << " < " << this->_abs_residual_tol);
            ML_CLOG(INFO, this->get_logger_id(), " Auk-b    Sweeper has converged w.r.t. absolute residual tolerance: "   << this->status()->abs_res_norm() << " < " << this->_abs_residual_tol);
          } else if (new_rel_res_norm_max < this->_rel_residual_tol) {
            ML_CLOG(INFO, this->get_logger_id(), "Sweeper has converged w.r.t. relative residual tolerance: "
                                               << new_rel_res_norm_max << " < " << this->_rel_residual_tol);
          } else {
            ML_CLOG(INFO, this->get_logger_id(), " uk+1- uk Sweeper has not yet converged to neither residual tolerance. "<< new_abs_res_norm_max << " < " << this->_abs_residual_tol);
            ML_CLOG(INFO, this->get_logger_id(), " Auk-b    Sweeper has not yet converged to neither residual tolerance. "<< this->status()->abs_res_norm() << " < " << this->_abs_residual_tol);
          }

          return (   new_abs_res_norm_max < this->_abs_residual_tol
                     || new_rel_res_norm_max < this->_rel_residual_tol);
        } else {
          ML_CLOG(WARNING, this->get_logger_id(), "No residual tolerances set. Thus skipping convergence check.");
          return false;
        }


        //////////////////////
      }
    }*/



  template<class SweeperTrait, class BaseFunction, typename Enabled>
  bool
  Sweeper<SweeperTrait, BaseFunction, Enabled>::converged()
  {
    return this->converged(false);
  }

  /**
   * @throws std::runtime_error in case the right node (i.e. the time end point) is not a quadrature node
   */
  template<class SweeperTrait, class BaseFunction,  typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::integrate_end_state(const typename SweeperTrait::time_t& dt)
  {
    UNUSED(dt);
    assert(this->get_quadrature() != nullptr);
    ML_CVLOG(4, this->get_logger_id(), "integrating end state");

    if (this->get_quadrature()->right_is_node()) {
      assert(this->get_end_state() != nullptr);
      assert(this->get_states().size() > 0);

      this->end_state()->data() = this->get_states().back()->get_data();
//       ML_CVLOG(1, this->get_logger_id(), "end state: " << to_string(this->get_end_state()));
    } else {
      throw std::runtime_error("integration of end state for quadrature not including right time interval boundary");
    }
  }

  /**
   * @throws std::runtime_error if not overwritten in specialized implementation
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::compute_residuals(const bool& only_last)
  {
    UNUSED(only_last);
    throw std::runtime_error("computation of residuals");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  Sweeper<SweeperTrait, BaseFunction, Enabled>::compute_residuals()
  {
    this->compute_residuals(false);
  }
}  // ::pfasst
