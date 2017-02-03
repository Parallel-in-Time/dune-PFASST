#include "pfasst/controller/controller.hpp"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
using std::shared_ptr;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT, class CommT>
  Controller<TransferT, CommT>::Controller()
    :   _status(std::make_shared<Status<typename TransferT::traits::fine_time_t>>())
      , _ready(false)
      , _logger_id("CONTROL")
  {}

  template<class TransferT, class CommT>
  shared_ptr<CommT>&
  Controller<TransferT, CommT>::communicator()
  {
    return this->_comm;
  }

  template<class TransferT, class CommT>
  const shared_ptr<CommT>
  Controller<TransferT, CommT>::get_communicator() const
  {
    return this->_comm;
  }

  template<class TransferT, class CommT>
  shared_ptr<Status<typename TransferT::traits::fine_time_t>>&
  Controller<TransferT, CommT>::status()
  {
    return this->_status;
  }

  template<class TransferT, class CommT>
  const shared_ptr<Status<typename TransferT::traits::fine_time_t>>
  Controller<TransferT, CommT>::get_status() const
  {
    return this->_status;
  }

  template<class TransferT, class CommT>
  size_t
  Controller<TransferT, CommT>::get_num_levels() const
  {
    return 0;
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::compute_num_steps()
  {
    if (this->get_status()->get_num_steps() != 0) {
      ML_CLOG(WARNING, this->get_logger_id(), "Total number of steps was already computed. Skipping.");
      return;
    }

    if (this->get_status()->get_t_end() <= 0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Time end point (" << this->get_status()->get_t_end()
                                         << ") must be non-zero positive.");
      throw std::logic_error("time end point must be non-zero positive");
    }

    if (this->get_status()->get_dt() <= 0.0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Time delta (" << this->get_status()->get_dt()
                                         << ") must be non-zero positive.");
      throw std::logic_error("time delta must be non-zero positive");
    }

    if (this->get_status()->get_time() >= this->get_status()->get_t_end()) {
      ML_CLOG(ERROR, this->get_logger_id(), "Time end point (" << this->get_status()->get_t_end()
                                         << ") must be greater than the current time point ("
                                         << this->get_status()->get_time() << ").");
      throw std::logic_error("time end point must be greater start time point");
    }

    const auto num_steps = (this->get_status()->get_t_end() - this->get_status()->get_time()) / this->get_status()->get_dt();
    ML_CLOG_IF(!almost_equal(num_steps * this->get_status()->get_dt(), lrint(num_steps) * this->get_status()->get_dt()),
            WARNING, this->get_logger_id(),
      LOG_FIXED << "End time point not an integral multiple of time delta: "
      << "(" << this->get_status()->get_t_end() << " - " << this->get_status()->get_time() << ") / " << this->get_status()->get_dt()
      << " = " << num_steps << " != " << lrint(num_steps));

    this->status()->num_steps() = lrint(num_steps);
  }

  template<class TransferT, class CommT>
  bool&
  Controller<TransferT, CommT>::ready()
  {
    return this->_ready;
  }

  template<class TransferT, class CommT>
  bool
  Controller<TransferT, CommT>::is_ready() const
  {
    return this->_ready;
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::set_logger_id(const std::string& logger_id)
  {
    this->_logger_id = logger_id;
  }

  template<class TransferT, class CommT>
  const char*
  Controller<TransferT, CommT>::get_logger_id() const
  {
    return this->_logger_id.c_str();
  }

  /**
   * @note Sets the maximum number of iterations and time end point from the command line arguments
   *   or leaves set values unchanged if not given on the command line.
   */
  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::set_options()
  {
    this->status()->max_iterations() = config::get_value<size_t>("num_iters", this->get_status()->get_max_iterations());
    this->status()->t_end() = config::get_value<typename TransferT::traits::fine_time_t>("t_end", this->get_status()->get_t_end());
  }

  template<class TransferT, class CommT>
  template<class SweeperT>
  void
  Controller<TransferT, CommT>::add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse)
  {
    UNUSED(sweeper); UNUSED(as_coarse);
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::add_transfer(shared_ptr<TransferT> transfer)
  {
    this->_transfer = transfer;
  }

  template<class TransferT, class CommT>
  const shared_ptr<TransferT>
  Controller<TransferT, CommT>::get_transfer() const
  {
    return this->_transfer;
  }

  template<class TransferT, class CommT>
  shared_ptr<TransferT>
  Controller<TransferT, CommT>::get_transfer()
  {
    return this->_transfer;
  }

  /**
   * @throws std::logic_error  if configured @p t_end is zero or negative.
   * @throws std::logic_error  if computed total number of steps times the configured time step
   *                           width @p dt does not lead to the desired @p t_end.
   */
  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::setup()
  {
    ML_CLOG_IF(this->is_ready(), WARNING, this->get_logger_id(),
      "Controller has already been setup.");

    ML_CVLOG(1, this->get_logger_id(), "setting up controller");

    if (this->get_status()->get_t_end() <= 0.0) {
      ML_CLOG(ERROR, this->get_logger_id(), "End time point must be larger than zero."
        << " (" << this->get_status()->get_t_end() << ")");
      throw std::logic_error("end time point must be larger zero");
    }

    this->compute_num_steps();
    const auto num_steps = this->get_status()->get_num_steps();
    if (!almost_equal(this->get_status()->get_time() + num_steps * this->get_status()->get_dt(),
                      this->get_status()->get_t_end())) {
      ML_CLOG(ERROR, this->get_logger_id(), "End time point not an integral multiple of time delta. " << LOG_FIXED
        << " (" << num_steps << " * " << this->get_status()->get_dt()
        << " = " << num_steps * this->get_status()->get_dt() << " != " << this->get_status()->get_t_end() << ")");
      throw std::logic_error("time end point is not an integral multiple of time delta");
    }

    ML_CLOG_IF(this->get_status()->get_max_iterations() == 0, WARNING, this->get_logger_id(),
      "You sould define a maximum number of iterations to avoid endless runs."
      << " (" << this->get_status()->get_max_iterations() << ")");

    this->ready() = true;
  }

  /**
   * @throws std::logic_error if `is_ready()` returns `false`
   */
  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::run()
  {
    if (!this->is_ready()) {
      ML_CLOG(ERROR, this->get_logger_id(), "Controller is not ready to run. setup() not called yet.");
      throw std::logic_error("controller not ready to run");
    }
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::post_run() {
    ML_CLOG(INFO, this->get_logger_id(), "Run Finished.");
  }

  template<class TransferT, class CommT>
  bool
  Controller<TransferT, CommT>::advance_time(const size_t& num_steps)
  {
    const time_t delta_time = num_steps * this->get_status()->get_dt();
    const time_t new_time = this->get_status()->get_time() + delta_time;

    auto status_summary = this->get_status()->summary();
    for (const auto& line : status_summary) {
      ML_CLOG(INFO, this->get_logger_id(), line);
    }

    if (new_time > this->get_status()->get_t_end() && !almost_equal(new_time, this->get_status()->get_t_end())) {
      ML_CLOG(WARNING, this->get_logger_id(), "Not advancing " << num_steps
                                           << ((num_steps > 1) ? " time steps " : " time step ")
                                           << "with dt=" << this->get_status()->get_dt() << " to t=" << new_time
                                           << " as it will exceed T_end=" << this->get_status()->get_t_end() << " by "
                                           << (new_time - this->get_status()->get_t_end()));

      return false;

    } else if(almost_equal(new_time, this->get_status()->get_t_end())) {
      ML_CLOG(INFO, this->get_logger_id(), "End time point reached: " << this->get_status()->get_t_end());

      return false;

    } else {
      ML_CVLOG(1, this->get_logger_id(), "Advancing " << num_steps
                                      << ((num_steps > 1) ? " time steps " : " time step ")
                                      << "with dt=" << this->get_status()->get_dt() << " to t=" << new_time);

      this->status()->time() += delta_time;
      this->status()->step() += num_steps;
      this->status()->iteration() = 0;

      return true;
    }
  }

  template<class TransferT, class CommT>
  bool
  Controller<TransferT, CommT>::advance_time()
  {
    return this->advance_time(1);
  }

  template<class TransferT, class CommT>
  bool
  Controller<TransferT, CommT>::advance_iteration()
  {
    if (this->get_status()->get_iteration() + 1 > this->get_status()->get_max_iterations()) {
      ML_CLOG(WARNING, this->get_logger_id(), "Not advancing to next iteration ("
                                           << (this->get_status()->get_iteration() + 1)
                                           << ") as it will exceed maximum number of allowed iterations ("
                                           << this->get_status()->get_max_iterations() << ")");

      return false;

    } else {
      ML_CVLOG(1, this->get_logger_id(), "Advancing to next iteration -> " << (this->get_status()->get_iteration() + 1));

      this->status()->iteration()++;

      return true;
    }
  }
}  // ::pfasst
