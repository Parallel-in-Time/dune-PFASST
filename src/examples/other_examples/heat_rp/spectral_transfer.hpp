#ifndef _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_
#define _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_

#include "pfasst/transfer/polynomial.hpp"

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/quadrature.hpp"
#include "pfasst/contrib/fft.hpp"


namespace pfasst
{
  namespace contrib
  {
    /**
     * @ingroup Contributed
     */
    template<
      class TransferTraits,
      typename Enabled = void
    >
    class SpectralTransfer
      : public PolynomialTransfer<TransferTraits>
    {
      public:
        using traits = TransferTraits;

      protected:
        pfasst::contrib::FFT<typename traits::fine_encap_t> fft;

      public:
        SpectralTransfer() = default;
        SpectralTransfer(const SpectralTransfer<TransferTraits, Enabled> &other) = default;
        SpectralTransfer(SpectralTransfer<TransferTraits, Enabled> &&other) = default;
        virtual ~SpectralTransfer() = default;
        SpectralTransfer<TransferTraits, Enabled>& operator=(const SpectralTransfer<TransferTraits, Enabled> &other) = default;
        SpectralTransfer<TransferTraits, Enabled>& operator=(SpectralTransfer<TransferTraits, Enabled> &&other) = default;

        virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                      shared_ptr<typename TransferTraits::fine_encap_t> fine);

        virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);
	
	virtual void restrict_u(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);
    };
  }  // ::pfasst::contrib
}  // ::pfasst

#include "sp.hpp"




#endif // _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_