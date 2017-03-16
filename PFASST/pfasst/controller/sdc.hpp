#ifndef _PFASST__CONTROLLER__SDC_HPP_
#define _PFASST__CONTROLLER__SDC_HPP_

#include <memory>
#include <type_traits>
using std::shared_ptr;

#include "pfasst/globals.hpp"
#include "pfasst/controller/status.hpp"
#include "pfasst/controller/controller.hpp"


namespace pfasst
{
  /**
   * @ingroup Controllers
   */
  template<
    class TransferT
  >
  class SDC
    : public Controller<TransferT>
  {
    static_assert(std::is_same<
                    std::integral_constant<size_t, TransferT::traits::num_levels>,
                    std::integral_constant<size_t, 1>>::value,
                  "SDC only works for single sweeper setups.");

    public:
      using transfer_t = TransferT;
      using time_t = typename transfer_t::traits::fine_time_t;

      static void init_loggers();

    protected:
      shared_ptr<typename transfer_t::traits::fine_sweeper_t> _sweeper;

      shared_ptr<void>               _transfer;
      shared_ptr<Status<time_t>>     _status;
      shared_ptr<comm::Communicator> _comm;
      bool                           _ready;

    public:
      SDC();
      SDC(const SDC<TransferT>& other) = default;
      SDC(SDC<TransferT>&& other) = default;
      virtual ~SDC() = default;
      SDC<TransferT>& operator=(const SDC<TransferT>& other) = default;
      SDC<TransferT>& operator=(SDC<TransferT>&& other) = default;

      virtual size_t get_num_levels() const override;

      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse);
      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper);

      virtual void add_transfer(shared_ptr<TransferT> transfer) override;

      virtual const shared_ptr<typename TransferT::traits::fine_sweeper_t> get_sweeper() const;
      virtual       shared_ptr<typename TransferT::traits::fine_sweeper_t> get_sweeper();

      virtual void set_options() override;

      virtual void setup() override;
      virtual void run() override;

      virtual bool advance_time(const size_t& num_steps) override;
      virtual bool advance_time() override;
      virtual bool advance_iteration() override;
  };
}  // ::pfasst

#include "pfasst/controller/sdc_impl.hpp"

#endif  // _PFASST__CONTROLLER__INTERFACE_HPP_
