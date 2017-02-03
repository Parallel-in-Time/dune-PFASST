#include "pfasst/controller/two_level_pfasst.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT, class CommT>
  TwoLevelPfasst<TransferT, CommT>::TwoLevelPfasst()
    : TwoLevelMLSDC<TransferT, CommT>()
  {
    TwoLevelPfasst<TransferT, CommT>::init_loggers();
    this->set_logger_id("PFASST");
    this->_prev_status = std::make_shared<Status<time_t>>();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::init_loggers()
  {
    log::add_custom_logger("PFASST");
    log::add_custom_logger("LVL_COARSE");
    log::add_custom_logger("LVL_FINE");
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::set_options()
  {
    TwoLevelMLSDC<TransferT, CommT>::set_options();
  }


  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::setup()
  {
    assert(this->get_communicator() != nullptr);

    TwoLevelMLSDC<TransferT, CommT>::setup();

    assert(this->get_transfer() != nullptr);

    if (this->get_num_levels() != 2) {
      ML_CLOG(ERROR, this->get_logger_id(), "Two levels (Sweeper) must have been added for Two-Level-PFASST.");
      throw std::logic_error("Two-Level-PFASST requires two levels");
    }

    if (this->get_communicator()->get_size() < 2) {
      ML_CLOG(ERROR, this->get_logger_id(), "Two-Level-PFASST requires at least two processes.");
      throw std::logic_error("two processes required for Two-Level-PFASST");
    }

    this->_prev_status = std::make_shared<Status<time_t>>();
    this->_prev_status->clear();
    this->_prev_status_temp = std::make_shared<Status<time_t>>();
    this->_prev_status_temp->clear();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::run()
  {
    Controller<TransferT, CommT>::run();

    assert(this->get_communicator() != nullptr);

    if (this->get_status()->get_num_steps() % this->get_communicator()->get_size() != 0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Number of time steps (" << this->get_status()->get_num_steps()
                                         << ") must be a multiple of the number of processors ("
                                         << this->get_communicator()->get_size() << ").");
      throw std::logic_error("number time steps must be multiple of number processors");
    }

    const size_t num_blocks = this->get_status()->get_num_steps() / this->get_communicator()->get_size();

    if (num_blocks == 0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Invalid Duration: There are more time processes ("
                                         << this->get_communicator()->get_size() << ") than time steps ("
                                         << this->get_status()->get_num_steps() << ").");
      throw std::logic_error("invalid duration: too many time processes for given time steps");
    }

    // iterate over time blocks (i.e. time-parallel blocks)
    do {
      this->status()->step() = this->_time_block * this->get_communicator()->get_size() \
                               + this->get_communicator()->get_rank();
      if (this->_time_block == 0) {
        this->status()->time() += this->get_status()->get_dt() * this->status()->get_step();
      }

      ML_CLOG(INFO, this->get_logger_id(), "");
      ML_CLOG(INFO, this->get_logger_id(), "Time Step " << (this->get_status()->get_step() + 1)
                                                        << " of " << this->get_status()->get_num_steps()
                                                        << " (i.e. t0=" << this->get_status()->get_time() << ")");

      // XXX: required here?
      this->_prev_status->clear();
      this->_prev_status_temp->clear();

      this->status()->set_primary_state(PrimaryState::PREDICTING);

      // iterate on each time step (i.e. iterations on single time step)
      do {
        if (this->get_status()->get_primary_state() == (+PrimaryState::PREDICTING)) {
          this->predictor();

        } else if (this->get_status()->get_primary_state() == (+PrimaryState::ITERATING)) {
          ML_CLOG(INFO, this->get_logger_id(), "");
          ML_CLOG(INFO, this->get_logger_id(), "Iteration " << this->get_status()->get_iteration());

          this->cycle_down();

          this->recv_coarse();
          this->sweep_coarse();
          this->send_coarse();

          this->cycle_up();

          this->sweep_fine();
          this->send_fine();

        } else {
          ML_CLOG(FATAL, this->get_logger_id(), "Something went severly wrong with the states.");
          ML_CLOG(FATAL, this->get_logger_id(), "Expected state: PREDICTING or ITERATING, got: "
                                                << (+this->get_status()->get_primary_state())._to_string());
          throw std::runtime_error("something went severly wrong");
        }

        // convergence check
      } while(this->advance_iteration());

      ML_CLOG(INFO, this->get_logger_id(), "");
      ML_CLOG(INFO, this->get_logger_id(), "Time Step done.");
    } while(this->advance_time(this->get_communicator()->get_size()));
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelPfasst<TransferT, CommT>::advance_time(const size_t& num_steps)
  {
    // receive potentially pending fine data of previous process
    this->recv_fine(true);
    this->get_communicator()->cleanup();

    if (TwoLevelMLSDC<TransferT, CommT>::advance_time(num_steps)) {
      ML_CLOG(INFO, this->get_logger_id(), "");
      this->broadcast();

      this->_time_block++;
      return true;
    } else {
      ML_CLOG(INFO, this->get_logger_id(), "");
      return false;
    }
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelPfasst<TransferT, CommT>::advance_time()
  {
    return this->advance_time(1);
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelPfasst<TransferT, CommT>::advance_iteration()
  {
    this->status()->set_primary_state(PrimaryState::INTER_ITER);

    this->get_check_prev_status();

    const bool fine_converged = this->get_fine()->converged(true);
    const bool previous_done = (this->get_communicator()->is_first())
                               ? true
                               : this->_prev_status->get_primary_state() <= (+PrimaryState::FAILED);
    ML_CLOG(DEBUG, this->get_logger_id(), "this status: " << to_string(this->status()));
    ML_CLOG(DEBUG, this->get_logger_id(), "prev status: " << to_string(this->_prev_status));

    if (previous_done && fine_converged) {
      ML_CLOG(INFO, this->get_logger_id(), "FINE sweeper has converged as well as previous process.");

      // receive potentially pending fine data of previous process
      this->recv_fine();

      this->status()->set_primary_state(PrimaryState::CONVERGED);
      this->_prev_status->set_primary_state(PrimaryState::UNKNOWN_PRIMARY);

    } else {
      ML_CLOG_IF(previous_done && !fine_converged, INFO, this->get_logger_id(),
        "previous process has converged but FINE sweeper not yet.");

      if (Controller<TransferT, CommT>::advance_iteration()) {
        ML_CLOG(INFO, this->get_logger_id(), "FINE sweeper has not yet converged and additional iterations to do.");
        this->get_fine()->save();
        this->get_coarse()->save();
        this->status()->set_primary_state(PrimaryState::ITERATING);

      } else {
        ML_CLOG(WARNING, this->get_logger_id(), "FINE sweeper has not yet converged and iterations threshold reached.");

        // receive potentially pending fine data of previous process
        this->recv_fine(true);

        this->status()->set_primary_state(PrimaryState::FAILED);
      }
    }

    this->send_status();

    this->get_fine()->converged(false);

    return (this->get_status()->get_primary_state() > (+PrimaryState::FAILED));
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::send_status()
  {
    if (!this->get_communicator()->is_last()) {
      ML_CVLOG(1, this->get_logger_id(), "sending status: " << to_string(this->get_status()));
      this->get_status()->send(this->get_communicator(),
                               this->get_communicator()->get_rank() + 1,
                               this->compute_tag(TagType::STATUS, TagLevel::FINE), false);
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::get_check_prev_status()
  {
    if (!this->get_communicator()->is_first()) {
      if (this->_prev_status->get_primary_state() > (+PrimaryState::FAILED)) {
        ML_CLOG(DEBUG, this->get_logger_id(), "prev known status: " << to_string(this->_prev_status));
        this->_prev_status_temp->clear();

        ML_CVLOG(1, this->get_logger_id(), "looking for updated state of previous process");
        this->_prev_status_temp->recv(this->get_communicator(),
                                      this->get_communicator()->get_rank() - 1,
                                      this->compute_tag(TagType::STATUS, TagLevel::FINE, TagModifier::PREV_STEP), true);
        // copy latest received status to the place where we use it from
        *(this->_prev_status) = *(this->_prev_status_temp);
        ML_CLOG(DEBUG, this->get_logger_id(), "Status received: " << to_string(this->_prev_status));

        if (this->_prev_status->get_primary_state() == (+PrimaryState::FAILED)) {
          ML_CLOG(WARNING, this->get_logger_id(), "previous process failed");

        } else if (this->_prev_status->get_primary_state() == (+PrimaryState::CONVERGED)) {
          ML_CLOG(WARNING, this->get_logger_id(), "previous process has converged; this process not");

        } else {
          ML_CVLOG(1, this->get_logger_id(), "previous process not finished");
        }
      }
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::recv_coarse()
  {
    if (!this->get_communicator()->is_first()) {
      if (this->_prev_status->get_primary_state() > (+PrimaryState::FAILED)) {
        ML_CVLOG(2, this->get_logger_id(), "looking for coarse data");
        this->get_coarse()
            ->initial_state()
            ->recv(this->get_communicator(),
                   this->get_communicator()->get_rank() - 1,
                   this->compute_tag(TagType::DATA, TagLevel::COARSE, TagModifier::PREV_STEP),
                   true);
      } else {
        ML_CLOG(WARNING, this->get_logger_id(), "previous process doesn't send any coarse data any more");
      }
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::send_coarse()
  {
    if (!this->get_communicator()->is_last()) {
      ML_CVLOG(1, this->get_logger_id(), "sending coarse end state");
      this->get_coarse()
          ->get_end_state()
          ->send(this->communicator(),
                 this->get_communicator()->get_rank() + 1,
                 this->compute_tag(TagType::DATA, TagLevel::COARSE),
                 true);
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::recv_fine(const bool& dummy)
  {
    if (!this->get_communicator()->is_first()) {
      ML_CVLOG(1, this->get_logger_id(), "looking for new initial value of fine level");
      const bool fine_avail = this->get_fine()
                                  ->initial_state()
                                  ->probe(this->get_communicator(),
                                          this->get_communicator()->get_rank() - 1,
                                          this->compute_tag(TagType::DATA, TagLevel::FINE,
                                                            (dummy)
                                                            ? TagModifier::PREV_STEP
                                                            : TagModifier::PREV_ITER_PREV_STEP));

      if (fine_avail) {
        this->get_fine()
            ->initial_state()
            ->recv(this->get_communicator(),
                   this->get_communicator()->get_rank() - 1,
                   this->compute_tag(TagType::DATA,
                                     TagLevel::FINE,
                                     (dummy)
                                     ? TagModifier::PREV_STEP
                                     : TagModifier::PREV_ITER_PREV_STEP),
                   true);
        ML_CVLOG(1, this->get_logger_id(), "new initial data on fine level received");
      } else {
        ML_CVLOG(1, this->get_logger_id(), "no new data available");
      }
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::send_fine()
  {
    if (!this->get_communicator()->is_last()) {
      ML_CVLOG(2, this->get_logger_id(), "sending fine data");
      this->get_fine()
          ->get_end_state()
          ->send(this->get_communicator(),
                 this->get_communicator()->get_rank() + 1,
                 this->compute_tag(TagType::DATA, TagLevel::FINE),
                 false);
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::cycle_down()
  {
    ML_CVLOG(1, this->get_logger_id(), "cycle down to coarse level");

    this->status()->set_secondary_state(SecondaryState::CYCLE_DOWN);

    this->get_transfer()->restrict(this->get_fine(), this->get_coarse(), true);
    this->get_transfer()->fas(this->get_status()->get_dt(), this->get_fine(), this->get_coarse());
    this->get_coarse()->save();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::cycle_up()
  {
    ML_CVLOG(1, this->get_logger_id(), "cycle up to fine level");

    this->status()->set_secondary_state(SecondaryState::CYCLE_UP);

    this->get_transfer()->interpolate(this->get_coarse(), this->get_fine(), true);

    this->recv_fine();

    this->get_transfer()->interpolate_initial(this->get_coarse(), this->get_fine());
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::predictor()
  {
    assert(this->get_status()->get_iteration() == 0);

    ML_CLOG(INFO, this->get_logger_id(), "");
    ML_CLOG(INFO, this->get_logger_id(), "Iteration 0 (PFASST Prediction)");

    // restrict fine initial condition ...
    this->get_transfer()->restrict_initial(this->get_fine(), this->get_coarse());
    // ... and spread it to all nodes on the coarse level
    this->get_coarse()->spread();
    this->get_coarse()->save();

    // perform PFASST prediction sweeps on coarse level
    for (size_t predict_step = 0;
         predict_step <= this->get_communicator()->get_rank();
         ++predict_step) {
      // do the sweeper's prediction once ...
      if (predict_step == 0) {
        this->predict_coarse();
      } else {
        // and default sweeps for subsequent processes
        this->recv_coarse();
        this->sweep_coarse();
      }

      this->send_coarse();
    }

    // return to fine level
    ML_CVLOG(1, this->get_logger_id(), "cycle up onto fine level");
    this->get_transfer()->interpolate(this->get_coarse(), this->get_fine(), true);
    this->sweep_fine();

    this->send_fine();

    // finalize prediction step
    this->get_coarse()->save();
    this->get_fine()->save();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::broadcast()
  {
    this->get_fine()->get_end_state()->bcast(this->get_communicator(), this->get_communicator()->get_size() - 1);
  }

  template<class TransferT, class CommT>
  int
  TwoLevelPfasst<TransferT, CommT>::compute_tag(const TagType type, const TagLevel level, const TagModifier mod) const
  {
    int tag = (type == (+TagType::DATA)) ? 1 : 0;

    if (type == (+TagType::DATA)) {
      const size_t iter = this->get_status()->get_iteration()
                          - ((   mod == (+TagModifier::PREV_ITER)
                              || mod == (+TagModifier::PREV_ITER_PREV_STEP)) ? 1 : 0);
      tag += (iter + 1) * 10000;
    }

    const size_t step = this->get_status()->get_step()
                        - ((   mod == (+TagModifier::PREV_STEP)
                            || mod == (+TagModifier::PREV_ITER_PREV_STEP)) ? 1 : 0);
    tag += (step + 1) * 1000000;

    tag += ((+level)._to_integral() + 1) * 100;

    ML_CLOG(DEBUG, this->get_logger_id(),
            "computing tag for " << (+type)._to_string() << " communication "
            << "on " << (+level)._to_string() << " level "
            << (mod == (+TagModifier::UNMOD) ? "without modifier" : "with modifier ")
            << (mod != (+TagModifier::UNMOD) ? (+mod)._to_string() : "")
            << " --> " << tag);

    return tag;
  }
}  // ::pfasst
