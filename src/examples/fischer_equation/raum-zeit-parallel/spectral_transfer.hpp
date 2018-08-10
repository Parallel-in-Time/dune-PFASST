#ifndef _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_
#define _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_

#include "pfasst/transfer/polynomial.hpp"

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/quadrature.hpp"
#include "pfasst/contrib/fft.hpp"

#include "../../datatypes/dune_vec.hpp"

const int dim=1;


namespace pfasst
{
  namespace contrib
  {
    /**
     * @ingroup Contributed
     */
    template<
      class TransferTraits
    >
    class SpectralTransfer
      : public PolynomialTransfer<TransferTraits>
    {
      public:
        using traits = TransferTraits;

      protected:
        pfasst::contrib::FFT<typename traits::fine_encap_t> fft;
        
        typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType;
        Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate_matrix;
        Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict_matrix;

      public:
        SpectralTransfer() = default;
        SpectralTransfer(const SpectralTransfer<TransferTraits> &other) = default;
        SpectralTransfer(SpectralTransfer<TransferTraits> &&other) = default;
        virtual ~SpectralTransfer() = default;
        SpectralTransfer<TransferTraits>& operator=(const SpectralTransfer<TransferTraits> &other) = default;
        SpectralTransfer<TransferTraits>& operator=(SpectralTransfer<TransferTraits> &&other) = default;

        /*virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                      shared_ptr<typename TransferTraits::fine_encap_t> fine);

        virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);
	
        virtual void restrict_u(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);*/
    
    
        virtual void set_matrix(Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate, Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict);

        virtual void create(std::shared_ptr<fe_manager> FinEl);
	
        virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                      shared_ptr<typename TransferTraits::fine_encap_t> fine);

        virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);

        virtual void restrict_u(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);
    
    
    
    };
  }  // ::pfasst::contrib
}  // ::pfasst







#include "spectral_transfer_impl.hpp"





#endif // _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_
