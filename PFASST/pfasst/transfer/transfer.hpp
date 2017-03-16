#ifndef _PFASST__TRANSFER__INTERFACE_HPP_
#define _PFASST__TRANSFER__INTERFACE_HPP_

#include <memory>
#include <type_traits>
using std::shared_ptr;

#include "pfasst/transfer/traits.hpp"


namespace pfasst
{
  /**
   * @defgroup TransferOp Transfer Operators
   *   Transfer Operators define the transition from one _level_ to another.
   * @ingroup Assistance
   */

  /**
   * @ingroup TransferOp
   */
  template<
    class TransferTraits,
    typename Enabled = void
  >
  class Transfer
    : public std::enable_shared_from_this<Transfer<TransferTraits, Enabled>>
  {
    public:
      using traits = TransferTraits;

      static_assert(std::is_convertible<
                      typename traits::coarse_time_t,
                      typename traits::fine_time_t
                    >::value,
                    "Coarse Time Type must be convertible to Fine Time Type");
      static_assert(std::is_convertible<
                      typename traits::fine_time_t,
                      typename traits::coarse_time_t
                    >::value,
                    "Fine Time Type must be convertible to Coarse Time Type");

    public:
      Transfer() = default;
      Transfer(const Transfer<TransferTraits, Enabled>& other) = default;
      Transfer(Transfer<TransferTraits, Enabled>&& other) = default;
      virtual ~Transfer() = default;
      Transfer<TransferTraits, Enabled>& operator=(const Transfer<TransferTraits, Enabled>& other) = default;
      Transfer<TransferTraits, Enabled>& operator=(Transfer<TransferTraits, Enabled>&& other) = default;

      virtual void interpolate_initial(const shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                                       shared_ptr<typename TransferTraits::fine_sweeper_t> fine);
      virtual void interpolate(const shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                               shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                               const bool initial = false);
      virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                    shared_ptr<typename TransferTraits::fine_encap_t> fine);

      virtual void restrict_initial(const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                    shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse);
      virtual void restrict(const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                            shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                            const bool initial = false);
      virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                 shared_ptr<typename TransferTraits::coarse_encap_t> coarse);

      virtual void fas(const typename TransferTraits::fine_time_t& dt,
                       const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                       shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse);
  };
}  // ::pfasst

#include "pfasst/transfer/transfer_impl.hpp"

#endif  // _PFASST__TRANSFER__INTERFACE_HPP_
