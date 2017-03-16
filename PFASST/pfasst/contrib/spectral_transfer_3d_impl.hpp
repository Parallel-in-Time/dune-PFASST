#include "pfasst/contrib/spectral_transfer.hpp"

#include <cassert>
#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/util.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"


namespace pfasst
{
  namespace contrib
  {
    /**
     * @internal
     * The interpolation is done in the Fourier space based on the data representation of FFTW3.
     *
     * The positive frequencies are located in the top-left corner of the 2D matrix, while the
     * negative frequencies are in the bottom-right corner:
     *
     * @verbatim
     * + + . .
     * + + . .
     * . . - -
     * . . - -
     * @endverbatim
     *
     * Interpolation is then simply inserting a "plus" of zeros in the center of the matrix of
     * frequencies keeping the general order of the frequencies:
     *
     * @verbatim
     * + + 0 0 0 0 . .
     * + + 0 0 0 0 . .
     * 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0
     * . . 0 0 0 0 - -
     * . . 0 0 0 0 - -
     * @endverbatim
     * @endinterl
     */
    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits,
      typename enable_if<
                 is_same<
                   typename TransferTraits::fine_encap_traits::dim_t,
                   integral_constant<size_t, 3>
                 >::value
               >::type>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                          shared_ptr<typename TransferTraits::fine_encap_t> fine)
    {
      ML_CVLOG(1, "TRANS", "interpolate data");

      const size_t coarse_ndofs = coarse->get_data().size();
      const size_t fine_ndofs = fine->get_data().size();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        fine->data() = coarse->get_data();

      } else {

      }
    }

    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits,
      typename enable_if<
                 is_same<
                   typename TransferTraits::fine_encap_traits::dim_t,
                   integral_constant<size_t, 3>
                 >::value
               >::type>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
    {
      ML_CVLOG(1, "TRANS", "restrict data");

      const size_t coarse_ndofs = coarse->get_data().size();
      const size_t fine_ndofs = fine->get_data().size();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        coarse->data() = fine->get_data();

      } else {
        const size_t coarse_dim_dofs = cbrt(coarse_ndofs);
        if (pow(coarse_dim_dofs, 3) != coarse_ndofs) {
          ML_CLOG(FATAL, "TRANS",
                  "Coarse space is not a cube: " << coarse_dim_dofs << "^3 != " << coarse_ndofs);
          throw runtime_error("coarse space not a cube");
        }
        const size_t fine_dim_dofs   = cbrt(fine_ndofs);
        if (pow(fine_dim_dofs, 3) != fine_ndofs) {
          ML_CLOG(FATAL, "TRANS",
                  "Fine space is not a cube: " << fine_dim_dofs << "^3 != " << fine_ndofs);
          throw runtime_error("fine space not a cube");
        }

        const size_t factor = fine_dim_dofs / coarse_dim_dofs;
        if (fine_dim_dofs != coarse_dim_dofs * 2) {
          ML_CLOG(FATAL, "TRANS", "FFTW based interpolation in 3D only for coarsening factor of 2");
          throw runtime_error("unsupported coarsening factor for FFTW interpolation");
        }

        for (size_t yi = 0; yi < coarse_dim_dofs; ++yi) {
          for (size_t xi = 0; xi < coarse_dim_dofs; ++xi) {
            for (size_t zi = 0; zi < coarse_dim_dofs; ++zi) {
              const size_t coarse_index = linearized_index(make_tuple(zi, yi, xi), coarse_dim_dofs);
              assert(coarse_index < coarse_ndofs);
              const size_t fine_index = factor * linearized_index(make_tuple(zi, yi, xi), fine_dim_dofs);
              assert(fine_index < fine_ndofs);
              coarse->data()[coarse_index] = fine->get_data()[fine_index];
            }
          }
        }
      }
    }
  }  // ::pfasst::contrib
}  // ::pfasst
