#include "pfasst/controller/sdc_n.hpp"

#include <memory>
#include <stdexcept>
#include <type_traits>
using std::shared_ptr;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT>
  SDC<TransferT>::SDC()
    : Controller<TransferT>()
  {
    SDC<TransferT>::init_loggers();
    this->set_logger_id("SDC");
  }

  template<class TransferT>
  void
  SDC<TransferT>::init_loggers()
  {
    log::add_custom_logger("SDC");
  }

  template<class TransferT>
  size_t
  SDC<TransferT>::get_num_levels() const
  {
    return ((this->_sweeper != nullptr) ? 1 : 0);
  }

  template<class TransferT>
  template<class SweeperT>
  void
  SDC<TransferT>::add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse)
  {
    UNUSED(as_coarse);
    this->add_sweeper(sweeper);
  }

  template<class TransferT>
  template<class SweeperT>
  void
  SDC<TransferT>::add_sweeper(shared_ptr<SweeperT> sweeper)
  {
    static_assert(std::is_same<SweeperT, typename TransferT::traits::fine_sweeper_t>::value,
                  "Sweeper must be a Fine Sweeper Type.");

    this->_sweeper = sweeper;
  }

  template<class TransferT>
  void
  SDC<TransferT>::add_transfer(shared_ptr<TransferT> transfer)
  {
    UNUSED(transfer);
    ML_CLOG(WARNING, this->get_logger_id(), "SDC Controller does not require a transfer operator.");
  }

  template<class TransferT>
  const shared_ptr<typename TransferT::traits::fine_sweeper_t>
  SDC<TransferT>::get_sweeper() const
  {
    return this->_sweeper;
  }

  template<class TransferT>
  shared_ptr<typename TransferT::traits::fine_sweeper_t>
  SDC<TransferT>::get_sweeper()
  {
    return this->_sweeper;
  }

  template<class TransferT>
  void
  SDC<TransferT>::set_options()
  {
    Controller<TransferT>::set_options();

    this->get_sweeper()->set_options();
  }

  template<class TransferT>
  void
  SDC<TransferT>::setup()
  {
    Controller<TransferT>::setup();

    if (this->get_num_levels() != 1) {
      ML_CLOG(ERROR, this->get_logger_id(), "One level (Sweeper) must have been added for SDC.");
      throw std::logic_error("SDC requires one level");
    }

    //std::cout <<  this->get_status().get_data().N() << std::endl; this->get_states().size()	
    this->get_sweeper()->status() = this->get_status();
    this->get_sweeper()->setup();
  }

  template<class TransferT>
  void
  SDC<TransferT>::run()
  {
    Controller<TransferT>::run();

    ML_CLOG(INFO, this->get_logger_id(), "");
    ML_CLOG(INFO, this->get_logger_id(), "Sequential SDC");
    ML_CLOG(INFO, this->get_logger_id(), "  t0:        " <<  this->get_status()->get_time());
    ML_CLOG(INFO, this->get_logger_id(), "  dt:        " <<  this->get_status()->get_dt());
    ML_CLOG(INFO, this->get_logger_id(), "  T:         " <<  this->get_status()->get_t_end());
    ML_CLOG(INFO, this->get_logger_id(), "  num steps: " <<  this->get_status()->get_num_steps());
    ML_CLOG(INFO, this->get_logger_id(), "  max iter:  " <<  this->get_status()->get_max_iterations());
//     ML_CLOG(INFO, this->get_logger_id(), "  Initial Value: " << to_string(this->get_sweeper()->get_initial_state()));

    // iterate over time steps
    do {
      ML_CLOG(INFO, this->get_logger_id(), "");
      ML_CLOG(INFO, this->get_logger_id(), "Time Step " << (this->get_status()->get_step() + 1)
                                        << " of " << this->get_status()->get_num_steps());

      // iterate on current time step
      do {
        const bool do_prediction = this->get_status()->get_iteration() == 0;

        if (do_prediction) {
	        
          ML_CLOG(INFO, this->get_logger_id(), "");
          ML_CLOG(INFO, this->get_logger_id(), "Iteration 0 (SDC Prediction)");
	  //std::cout << "-----------------------------------------------------------------------------------------                vorm pre_predict" <<  std::endl;
          this->get_sweeper()->pre_predict();

	  //std::cout << "-----------------------------------------------------------------------------------------                vorm predict" <<  std::endl;
	  this->get_sweeper()->predict();
	  	        
	  //std::cout << "-----------------------------------------------------------------------------------------                vorm post_predict" <<  std::endl;
          this->get_sweeper()->post_predict();

          this->get_sweeper()->pre_sweep();
	  //std::cout << "-----------------------------------------------------------------------------------------                vorm sweep" <<  std::endl;
          this->get_sweeper()->sweep();
	  //std::cout << "-----------------------------------------------------------------------------------------                vorm post_sweep" <<  std::endl;
          this->get_sweeper()->post_sweep();
	        

        } else {
	        
          ML_CLOG(INFO, this->get_logger_id(), "");
          ML_CLOG(INFO, this->get_logger_id(), "Iteration " << this->get_status()->get_iteration());
	  //std::cout << "-----------------------------------------------------------------------------------------                vorm presweep" <<  std::endl;
          this->get_sweeper()->pre_sweep();
	  //std::cout << "-----------------------------------------------------------------------------------------                vorm sweep" <<  std::endl;
          this->get_sweeper()->sweep();
	  //std::cout << "-----------------------------------------------------------------------------------------                vorm post_sweep" <<  std::endl;
          this->get_sweeper()->post_sweep();
        }
      } while(this->advance_iteration());
      
//       stringstream toss;
//       toss << this->get_status()->get_step();
//       string name = toss.str();
//       //if(this->get_sweeper()._write) 
//           this->get_sweeper()->write_results(this->get_sweeper()->get_end_state(), name); 
      
    } while(this->advance_time());
  }

  template<class TransferT>
  bool
  SDC<TransferT>::advance_time(const size_t& num_steps)
  {
    this->get_sweeper()->post_step();

    if (Controller<TransferT>::advance_time(num_steps)) {
      this->get_sweeper()->advance(num_steps);
      return true;
    } else {
      return false;
    }
  }

  template<class TransferT>
  bool
  SDC<TransferT>::advance_time()
  {
    return this->advance_time(1);
  }

  template<class TransferT>
  bool
  SDC<TransferT>::advance_iteration()
  {
    if (this->get_sweeper()->converged(false)) {
      ML_CLOG(INFO, this->get_logger_id(), "Sweeper has converged.");
      return false;
    } else if (Controller<TransferT>::advance_iteration()) {
      ML_CVLOG(1, this->get_logger_id(), "Sweeper has not yet converged and additional iterations to do.");
      this->get_sweeper()->save();
      return true;
    } else {
      ML_CLOG(WARNING, this->get_logger_id(), "Sweeper has not yet converged and no more iterations to do.");
      return false;
    }
  }
}  // ::pfasst
