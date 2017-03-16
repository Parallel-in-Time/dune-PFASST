#ifndef _PFASST__CONTROLLER__TWO_LEVEL_MLSDC_HPP_
#define _PFASST__CONTROLLER__TWO_LEVEL_MLSDC_HPP_

#include <memory>
using std::shared_ptr;

#include "pfasst/controller/controller.hpp"
#include "pfasst/comm/communicator.hpp"


namespace pfasst
{
  /**
   * @ingroup Controllers
   */
  template<
    class TransferT,
    class CommT = comm::Communicator
  >
  class TwoLevelMLSDC
    : public Controller<TransferT, CommT>
  {
    public:
      using transfer_t = TransferT;
      using comm_t = CommT;
      using time_t = typename transfer_t::traits::fine_time_t;

      static void init_loggers();

    protected:
      shared_ptr<typename transfer_t::traits::coarse_sweeper_t> _coarse_level;
      shared_ptr<typename transfer_t::traits::fine_sweeper_t>   _fine_level;

      virtual void predict_coarse();
      virtual void predict_fine();
      virtual void sweep_coarse();
      virtual void sweep_fine();

      virtual void cycle_down();
      virtual void cycle_up();

    public:
      TwoLevelMLSDC();
      TwoLevelMLSDC(const TwoLevelMLSDC<TransferT, CommT>& other) = default;
      TwoLevelMLSDC(TwoLevelMLSDC<TransferT, CommT>&& other) = default;
      virtual ~TwoLevelMLSDC() = default;
      TwoLevelMLSDC<TransferT, CommT>& operator=(const TwoLevelMLSDC<TransferT, CommT>& other) = default;
      TwoLevelMLSDC<TransferT, CommT>& operator=(TwoLevelMLSDC<TransferT, CommT>&& other) = default;

      virtual size_t get_num_levels() const override;

      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse);

      virtual const shared_ptr<typename TransferT::traits::coarse_sweeper_t> get_coarse() const;
      virtual       shared_ptr<typename TransferT::traits::coarse_sweeper_t> get_coarse();
      virtual const shared_ptr<typename TransferT::traits::fine_sweeper_t> get_fine() const;
      virtual       shared_ptr<typename TransferT::traits::fine_sweeper_t> get_fine();

      virtual void set_options() override;

      virtual void setup() override;
      virtual void run() override;

      virtual bool advance_time(const size_t& num_steps) override;
      virtual bool advance_time() override;
      virtual bool advance_iteration() override;
  };
}  // ::pfasst

#include "pfasst/controller/two_level_mlsdc_impl.hpp"

#endif  // _PFASST__CONTROLLER__TWO_LEVEL_MLSDC_HPP_
