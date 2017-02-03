#ifndef _PFASST__CONTROLLER__INTERFACE_HPP_
#define _PFASST__CONTROLLER__INTERFACE_HPP_

#include <memory>
#include <string>
using std::shared_ptr;

#include "pfasst/globals.hpp"
#include "pfasst/comm/communicator.hpp"
#include "pfasst/controller/status.hpp"


namespace pfasst
{
  /**
   * The Controller is the main driver to execute a Sweeper on multiple time steps.
   *
   * @tparam TransferT  Type of the transfer operator between two sweepers.
   * @tparam CommT      Type of the Communicator for parallel execution.
   *
   * @ingroup Controllers
   */
  template<
    class TransferT,
    class CommT = comm::Communicator
  >
  class Controller
    : public std::enable_shared_from_this<Controller<TransferT, CommT>>
  {
    public:
      //! Type of the Transfer operator.
      using transfer_t = TransferT;
      //! Type of the Communicator.
      using comm_t = CommT;
      //! Type of the time domain.
      using time_t = typename transfer_t::traits::fine_time_t;

    protected:
      //! Communicator instance.
      shared_ptr<comm_t>          _comm;

      //! Transfer operator instance.
      shared_ptr<transfer_t>      _transfer;
      //! Status object.
      shared_ptr<Status<time_t>>  _status;
      //! Flag to indicate readiness for execution.
      bool                        _ready;
      //! Name of the Controller in the logs.
      std::string                 _logger_id;

      /**
       * Compute total number of steps.
       *
       * Given the initial time point, the desired end time point and width of each time step,
       * compute the total number of time steps required to propagate the solution to the end time
       * point.
       */
      virtual void compute_num_steps();
      /**
       * Accessor for the readiness.
       *
       * @returns readiness flag
       */
      virtual bool& ready();

    public:
      //! @{
      Controller();
      Controller(const Controller<TransferT, CommT>& other) = default;
      Controller(Controller<TransferT, CommT>&& other) = default;
      virtual ~Controller() = default;
      Controller<TransferT, CommT>& operator=(const Controller<TransferT, CommT>& other) = default;
      Controller<TransferT, CommT>& operator=(Controller<TransferT, CommT>&& other) = default;
      //! @}

      //! @name Accessors
      //! @{
      /**
       * Accessor for the Communicator.
       */
      virtual       shared_ptr<CommT>& communicator();
      //! Read-only version of `communicator()`.
      virtual const shared_ptr<CommT>  get_communicator() const;

      /**
       * Accessor for the Status.
       */
      virtual       shared_ptr<Status<typename TransferT::traits::fine_time_t>>& status();
      //! Read-only version of `status()`.
      virtual const shared_ptr<Status<typename TransferT::traits::fine_time_t>>  get_status() const;

      /**
       * Number of levels/sweeper currently configured for this Controller.
       *
       * Each call to `Controller::add_sweeper` will increase this count by one.
       */
      virtual size_t get_num_levels() const;
      /**
       * Read-only accessor for the readiness of the Controller.
       *
       * @returns `true` when `setup()` was called successfully, `false` otherwise.
       */
      virtual bool   is_ready() const;

      /**
       * Setter for the logger ID.
       *
       * @param[in] logger_id  new name of the Controller instance in the logs.
       */
      virtual       void  set_logger_id(const std::string& logger_id);
      /**
       * Read-only accessor for the logger ID.
       */
      virtual const char* get_logger_id() const;

      /**
       * Add another Sweeper as a new level.
       *
       * @tparam SweeperT  Type of the Sweeper. Usually, this is inferred by the first argument.
       *
       * @param[in] sweeper    pointer to the Sweeper instance
       * @param[in] as_coarse  if `true`, the @p sweeper is added as the currently coarsest;
       *                       otherwise it's added as the currently finest Sweeper
       */
      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse);

      /**
       * Add/set the Transfer Operator.
       *
       * @param[in] transfer pointer to the Transfer Operator instance
       *
       * @todo Consider removing this from the base Controller in favour of a more cleaner interface.
       */
      virtual void add_transfer(shared_ptr<TransferT> transfer);

      /**
       * Accessor for the Transfer Operator.
       *
       * @todo Rename non-const `Controller::get_transfer()` to `Controller::transfer()` to match
       *    other non-const accessors.
       */
      virtual       shared_ptr<TransferT> get_transfer();
      //! Read-only version of `get_transfer()`.
      virtual const shared_ptr<TransferT> get_transfer() const;
      //! @}

      //! @name Configuration and Setup
      //! @{
      /**
       * Configures the Controller according to given command line arguments.
       */
      virtual void set_options();
      /**
       * Central Setup routine to prepare the Controller, Sweepers and Transfer Operators.
       *
       * This function has to be called exactly once before `run()`.
       */
      virtual void setup();
      //! @}

      //! @name Execution
      //! @{
      /**
       * Central algorithm to drive the Sweepers.
       *
       * In this base implementation it does nothing but checking the readiness of the Controller.
       * The actual logic has to be implemented in specializations.
       */
      virtual void run();
      /**
       * Post Run Hook.
       *
       * Might get overwritten by specializations to do some basic logic after the whole execution
       * is finished.
       */
      virtual void post_run();
      //! @}

      //! @name Utilities
      //! @{
      /**
       * Advances given number of steps.
       *
       * The current start time point is extended by the configured time step width.
       * The resulting time point is compared against the configured time end point.
       *
       * @param[in] num_steps  number of steps to advance
       *
       * @returns `true` if the next time step is required to reach desired time end point. `false`
       *   otherwise.
       */
      virtual bool advance_time(const size_t& num_steps);
      /**
       * Advance to the very next time step.
       *
       * @returns `advance_time(const size_t&)`.
       *
       * @note Calls `advance_time(const size_t&)` with `1`.
       *
       * @overload
       */
      virtual bool advance_time();
      /**
       * Advance to next iteration in the current time step.
       *
       * Checks the current iteration index w.r.t. to the total number of allowed iterations (c.f.
       * `Status::get_max_iterations()`) and determines whether an additional iteration is allowed.
       *
       * @returns `true` if the next iteration should be done, `false` otherwise
       */
      virtual bool advance_iteration();
      //! @}
  };
}  // ::pfasst

#include "pfasst/controller/controller_impl.hpp"

#endif  // _PFASST__CONTROLLER__INTERFACE_HPP_
