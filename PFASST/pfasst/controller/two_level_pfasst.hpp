#ifndef _PFASST__CONTROLLER__TWO_LEVEL_PFASST_HPP_
#define _PFASST__CONTROLLER__TWO_LEVEL_PFASST_HPP_

#include <memory>
using std::shared_ptr;


#include <better-enums/enum.h>


#include "pfasst/controller/two_level_mlsdc.hpp"
#include "pfasst/comm/mpi_p2p.hpp"


namespace pfasst
{
  namespace detail
  {
    ENUM(TagLevel, size_t, FINE, COARSE, ANY)
    ENUM(TagModifier, size_t, PREV_STEP, PREV_ITER, PREV_ITER_PREV_STEP, UNMOD)
    ENUM(TagType, size_t, STATUS, DATA)
  }  // ::pfasst::detail


  /**
   * @ingroup Controllers
   */
  template<
    class TransferT,
    class CommT = comm::MpiP2P
  >
  class TwoLevelPfasst
    : public TwoLevelMLSDC<TransferT, CommT>
  {
    using TagLevel = pfasst::detail::TagLevel;
    using TagModifier = pfasst::detail::TagModifier;
    using TagType = pfasst::detail::TagType;

    public:
      using transfer_t = TransferT;
      using comm_t = CommT;
      using time_t = typename transfer_t::traits::fine_time_t;

      static void init_loggers();

    protected:
      shared_ptr<Status<time_t>> _prev_status;
      shared_ptr<Status<time_t>> _prev_status_temp;
      size_t _time_block = 0;

      virtual void send_status();
      virtual void get_check_prev_status();

      virtual void recv_coarse();
      virtual void send_coarse();
      virtual void recv_fine(const bool& dummy = false);
      virtual void send_fine();

      virtual void predictor();
      virtual void cycle_down() override;
      virtual void cycle_up() override;

      virtual void broadcast();

      int compute_tag(const TagType type,
                      const TagLevel level = TagLevel::ANY,
                      const TagModifier mod = TagModifier::UNMOD) const;

    public:
      TwoLevelPfasst();
      TwoLevelPfasst(const TwoLevelPfasst<TransferT, CommT>& other) = default;
      TwoLevelPfasst(TwoLevelPfasst<TransferT, CommT>&& other) = default;
      virtual ~TwoLevelPfasst() = default;
      TwoLevelPfasst<TransferT, CommT>& operator=(const TwoLevelPfasst<TransferT, CommT>& other) = default;
      TwoLevelPfasst<TransferT, CommT>& operator=(TwoLevelPfasst<TransferT, CommT>&& other) = default;

      virtual void set_options() override;

      virtual void setup() override;
      virtual void run() override;

      virtual bool advance_time(const size_t& num_steps) override;
      virtual bool advance_time() override;
      virtual bool advance_iteration() override;
  };
}  // ::pfasst

#include "pfasst/controller/two_level_pfasst_impl.hpp"

#endif  // _PFASST__CONTROLLER__TWO_LEVEL_PFASST_HPP_
