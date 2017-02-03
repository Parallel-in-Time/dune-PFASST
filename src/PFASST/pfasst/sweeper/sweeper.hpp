#ifndef _PFASST__SWEEPER__INTERFACE_HPP_
#define _PFASST__SWEEPER__INTERFACE_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/sweeper/traits.hpp"
#include "pfasst/controller/status.hpp"
#include "pfasst/encap/encapsulation.hpp"
#include "pfasst/quadrature.hpp"
using pfasst::quadrature::IQuadrature;


namespace pfasst
{
  /**
   * @defgroup Sweepers Sweepers
   *   Sweepers represent the core method for advancing from one time point to another.
   * @ingroup Controllers
   */

  /**
   * Base Sweeper defining the generic interface for advancing to the next time step.
   *
   * @note As a general rule of thumb, whenever there is a `std::vector` of Encapsulation values, it
   *   is the spatial data at all time nodes including the initial time point.
   *   Thus, these vectors have one element more than the total number of quadrature nodes.
   *   This additional element is the very first one.
   *
   * @note It is assumed, that only Sweepers based on Spectral Deferred Corrections are implemented.
   *
   * @tparam SweeperTrait  Type Trait for the sweeper providing type information on used data types.
   *   The following is asserted at compile time:
   *     * `std::is_arithmetic<SweeperTrait::time_t>`
   *     * `std::is_constructible<SweeperTrait::encap_t>`
   *     * `std::is_destructible<SweeperTrait::encap_t>`
   * @tparam Enabled       utility parameter used for partial specialization
   *
   * @ingroup Sweepers
   */
  template<
    class SweeperTrait,
    typename Enabled = void
  >
  class Sweeper
    : public std::enable_shared_from_this<Sweeper<SweeperTrait, Enabled>>
  {
    public:
      /**
       * Sweeper type traits providing important type information for the sweeper.
       *
       * @see pfasst::sweeper_traits
       */
      using traits = SweeperTrait;

    static_assert(std::is_arithmetic<typename traits::time_t>::value,
                  "precision type must be arithmetic");
    static_assert(std::is_constructible<typename traits::encap_t>::value,
                  "Encapsulation type must be constructible");
    static_assert(std::is_destructible<typename traits::encap_t>::value,
                  "Encapsulation type must be destructible");

    protected:
      //! Quadrature used by the Sweeper.
      shared_ptr<IQuadrature<typename traits::time_t>>    _quadrature;
      //! Factory for instantiating Encapsulation objects for spatial data.
      shared_ptr<typename traits::encap_t::factory_t>     _factory;

      //! Vector of the spatial data at all nodes including the initial value.
      vector<shared_ptr<typename traits::encap_t>>        _states;
      //! Vector of the spatial data at all nodes of the previous iteration.
      vector<shared_ptr<typename traits::encap_t>>        _previous_states;
      //! Spatial data at the desired end time point.
      shared_ptr<typename traits::encap_t>                _end_state;

      //! Values for the coarse level correction at all time nodes.
      vector<shared_ptr<typename traits::encap_t>>        _tau;
      //! Full spatial residual at all time nodes.
      vector<shared_ptr<typename traits::encap_t>>        _residuals;
      //! Norms of absolute residuals at all time nodes.
      vector<typename traits::spatial_t>                  _abs_res_norms;
      //! Norms of relative (w.r.t. previous iteration) residuals at all time nodes.
      vector<typename traits::spatial_t>                  _rel_res_norms;

      //! Status object.
      shared_ptr<Status<typename traits::time_t>>         _status;
      //! Tolerance for the absolute residual.
      typename traits::spatial_t                          _abs_residual_tol;
      //! Tolerance for the relative residual.
      typename traits::spatial_t                          _rel_residual_tol;

      /**
       * Name of the Sweeper in the logs.
       *
       * Defaults to `SWEEPER`.
       */
      string                                              _logger_id;

      //! @{
      /**
       * Compute the solution at the desired end time point.
       *
       * In cases where the quadrature nodes include the desired end time point this does a simple
       * copy of the last state (`states().back()`) to the end state (`end_state()`).
       * Otherwise, an additional integration step is conducted.
       *
       * @param[in] dt don't know why this is there ...
       *
       * @todo Do we need this @p dt ?
       */
      virtual void integrate_end_state(const typename SweeperTrait::time_t& dt);
      /**
       * Compute full spatial residuals.
       *
       * @param[in] only_last if `true` computes only the residual at the desired end time point
       */
      virtual void compute_residuals(const bool& only_last);
      /**
       * Computes full spatial residuals at all time nodes.
       *
       * Calls `compute_residuals(const bool&)` with `false`.
       *
       * @overload
       */
      virtual void compute_residuals();
      //! @}

      //! @{
      /**
       * Accessor for the spatial data of the previous iteration at all time nodes.
       */
      virtual vector<shared_ptr<typename SweeperTrait::encap_t>>& previous_states();
      /**
       * Accessor for the spatial data at the end point.
       */
      virtual shared_ptr<typename SweeperTrait::encap_t>&         end_state();
      /**
       * Accessor for the full spatial residuals at all time nodes.
       */
      virtual vector<shared_ptr<typename SweeperTrait::encap_t>>& residuals();
      //! @}

      //! @{
      /**
       * Routine to instantiate all the internal variables.
       */
      virtual void initialize();
      //! @}

    public:
      //! @{
      explicit Sweeper();
      Sweeper(const Sweeper<SweeperTrait, Enabled>& other) = default;
      Sweeper(Sweeper<SweeperTrait, Enabled>&& other) = default;
      virtual ~Sweeper() = default;
      Sweeper<SweeperTrait, Enabled>& operator=(const Sweeper<SweeperTrait, Enabled>& other) = default;
      Sweeper<SweeperTrait, Enabled>& operator=(Sweeper<SweeperTrait, Enabled>&& other) = default;
      //! @}

      //! @name Accessors
      //! @{
      /**
       * Accessor to the used Quadrature rule.
       */
      virtual       shared_ptr<IQuadrature<typename SweeperTrait::time_t>>& quadrature();
      //! Read-only version of `quadrature()`.
      virtual const shared_ptr<IQuadrature<typename SweeperTrait::time_t>>  get_quadrature() const;

      /**
       * Accessor for the Status object.
       */
      virtual       shared_ptr<Status<typename SweeperTrait::time_t>>& status();
      //! Read-only version of `status()`.
      virtual const shared_ptr<Status<typename SweeperTrait::time_t>>  get_status() const;

      /**
       * Accessor for the EncapsulationFactory used to initialize new spatial data objects.
       */
      virtual       shared_ptr<typename SweeperTrait::encap_t::factory_t>& encap_factory();
      //! Read-only version of `encap_factory()`.
      virtual const typename SweeperTrait::encap_t::factory_t&             get_encap_factory() const;

      /**
       * Accessor for the spatial data at the initial time point.
       *
       * @note This is a shortcut for `get_states().front()` with additional error checking.
       */
      virtual       shared_ptr<typename SweeperTrait::encap_t>&         initial_state();
      //! Read-only version of `initial_state()`.
      virtual const shared_ptr<typename SweeperTrait::encap_t>          get_initial_state() const;

      /**
       * Accessor for the spatial data at all time points.
       */
      virtual       vector<shared_ptr<typename SweeperTrait::encap_t>>& states();
      //! Read-only version of `states()`.
      virtual const vector<shared_ptr<typename SweeperTrait::encap_t>>& get_states() const;

      /**
       * Accessor for the coarse level correction data at all time points.
       */
      virtual       vector<shared_ptr<typename SweeperTrait::encap_t>>& tau();
      //! Read-only version of `tau()`.
      virtual const vector<shared_ptr<typename SweeperTrait::encap_t>>& get_tau() const;

      /**
       * Accessor for the spatial data at all time points of the previous iteration.
       */
      virtual const vector<shared_ptr<typename SweeperTrait::encap_t>>& get_previous_states() const;

      /**
       * Read-only accessor for the spatial data at the time step's end point.
       *
       * @note This will get computed and set by integrate_end_state().
       */
      virtual const shared_ptr<typename SweeperTrait::encap_t>          get_end_state() const;

      /**
       * Read-only accessor for the residuals at each time point.
       */
      virtual const vector<shared_ptr<typename SweeperTrait::encap_t>>& get_residuals() const;

      /**
       * Setter for the logger ID.
       *
       * @param[in] logger_id  new name of the Sweeper instance in the logs.
       */
      virtual       void  set_logger_id(const string& logger_id);
      /**
       * Read-only accessor for the logger ID.
       */
      virtual const char* get_logger_id() const;
      //! @}

      //! @name Configuration and Setup
      //! @{
      /**
       * Configures the Sweeper according to given command line arguments.
       */
      virtual void set_options();
      /**
       * Sets the tolerance for the absolute residual to the given value.
       *
       * @note This tolerance is used by `converged()` as a termination criterion.
       *
       * @param[in] abs_res_tol new tolerance for the absolute residual
       */
      virtual void set_abs_residual_tol(const typename SweeperTrait::spatial_t& abs_res_tol);
      /**
       * Sets the tolerance for the relative residual to the given value.
       *
       * @note This tolerance is used by `converged()` as a termination criterion.
       *
       * @param[in] rel_res_tol new tolerance for the relative residual
       */
      virtual void set_rel_residual_tol(const typename SweeperTrait::spatial_t& rel_res_tol);

      /**
       * Generic setup routine for the Sweeper.
       *
       * This function does basic consistency checking on the configuration and calls `initialize()`
       * to prepare the Sweeper for the actual execution.
       */
      virtual void setup();
      /**
       * Resets all internal states and variables of the Sweeper.
       */
      virtual void reset();
      //! @}

      //! @name Prediction Step
      //! @{
      /**
       * Pre-Prediction Step.
       *
       * Might get overwritten by specializations to do some basic logic prior to the prediction
       * step.
       */
      virtual void pre_predict();
      /**
       * Prediction Step.
       *
       * In most cases, the very first iteration requires special handling and is called a
       * _prediction step_.
       */
      virtual void predict();
      /**
       * Post-Prediction Step.
       *
       * Might get overwritten by specializations to do some basic logic at the end of each
       * iteration.
       *
       * @note Calls `integrate_end_state()`.
       *   When overwriting, make sure to either call this overwritten function or call
       *   `integrate_end_state()`.
       */
      virtual void post_predict();
      //! @}

      //! @name Sweep Step
      //! @{
      /**
       * Pre-Sweep Step.
       *
       * Might get overwritten by specializations to do some basic logic prior to each sweep step.
       */
      virtual void pre_sweep();
      /**
       * Sweep Step.
       *
       * The central function implementing the sweeping algorithm to compute the solution at all
       * time nodes.
       */
      virtual void sweep();
      /**
       * Post-Sweep Step.
       *
       * Might get overwritten by specializations to do some basic logic at the end of each
       * iteration.
       *
       * @note Calls `integrate_end_state()`.
       *   When overwriting, make sure to either call this overwritten function or call
       *   `integrate_end_state()`.
       */
      virtual void post_sweep();
      //! @}

      //! @name Post-Step and Utilities
      //! @{
      /**
       * Post-Step.
       *
       * Might get overwritten by specializations to do some basic logic after each time step.
       */
      virtual void post_step();

      /**
       * Advance to the next time step.
       *
       * @param[in] num_steps advance given number of time steps.
       */
      virtual void advance(const size_t& num_steps);
      /**
       * Advance a to the very next time step.
       *
       * @note Calls `advance(const size_t&)` with `1`.
       *
       * @overload
       */
      virtual void advance();

      /**
       * Spread the initial value to all time nodes.
       */
      virtual void spread();
      /**
       * Save solution at all time nodes.
       *
       * This copies the spatial values of all `states()` to `previous_states()`.
       */
      virtual void save();
      /**
       * Evaluate the right hand side of the problem equation.
       *
       * @param[in] initial_only if `true` only at the initial value
       */
      virtual void reevaluate(const bool initial_only);
      /**
       * Evaluate the right hand side of the problem equation at all time nodes.
       *
       * @note calls `reevaluate(const bool)` with `false`.
       *
       * @overload
       */
      virtual void reevaluate();
      /**
       * Integrates spatial data at all time nodes for given time step width.
       *
       * @param[in] dt width of time interval.
       */
      virtual vector<shared_ptr<typename SweeperTrait::encap_t>> integrate(const typename SweeperTrait::time_t& dt);

      /**
       * Tests current solution residuals for tolerances.
       *
       * @param[in] pre_check if `true` only tests residuals at the last time point.
       *
       * @returns `true` when either maximum of the relative or absolute residual norm at all time
       *   nodes is below the configured tolerances. `false` otherwise.
       */
      virtual bool converged(const bool pre_check);
      /**
       * Tests residuals at the end time point against configured tolerances.
       *
       * @note calls `converged(const bool)` with `false`.
       *
       * @overload
       */
      virtual bool converged();
      //! @}
  };
}

#include "pfasst/sweeper/sweeper_impl.hpp"

#endif  // _PFASST__SWEEPER__INTERFACE_HPP_
