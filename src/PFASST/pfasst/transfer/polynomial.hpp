#ifndef _PFASST__TRANSFER__POLYNOMIAL_HPP_
#define _PFASST__TRANSFER__POLYNOMIAL_HPP_

#include "pfasst/transfer/transfer.hpp"

#include <memory>
using std::shared_ptr;

#include "pfasst/globals.hpp"
#include "pfasst/quadrature.hpp"


namespace pfasst
{
  /**
   * @ingroup TransferOp
   */
  template<
    class TransferTraits,
    typename Enabled = void
  >
  class PolynomialTransfer
    : public Transfer<TransferTraits, Enabled>
  {
    public:
      using traits = TransferTraits;

    protected:
      Matrix<typename traits::fine_time_t> tmat;
      Matrix<typename traits::fine_time_t> fmat;

      virtual void setup_tmat(const shared_ptr<quadrature::IQuadrature<typename TransferTraits::fine_time_t>> fine_quad,
                              const shared_ptr<quadrature::IQuadrature<typename TransferTraits::coarse_time_t>> coarse_quad);

    public:
      PolynomialTransfer() = default;
      PolynomialTransfer(const PolynomialTransfer<TransferTraits, Enabled>& other) = default;
      PolynomialTransfer(PolynomialTransfer<TransferTraits, Enabled>&& other) = default;
      virtual ~PolynomialTransfer() = default;
      PolynomialTransfer<TransferTraits, Enabled>& operator=(const PolynomialTransfer<TransferTraits, Enabled>& other) = default;
      PolynomialTransfer<TransferTraits, Enabled>& operator=(PolynomialTransfer<TransferTraits, Enabled>&& other) = default;

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

#include "pfasst/transfer/polynomial_impl.hpp"

#endif  // _PFASST__TRANSFER__POLYNOMIAL_HPP_
