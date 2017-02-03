#ifndef _PFASST__CONTRIB__FFT_HPP_
#define _PFASST__CONTRIB__FFT_HPP_

#include <complex>
using std::real;
#include <map>
#include <memory>
#include <utility>
using namespace std;


//#include <fftw3.h>



namespace pfasst
{
  /**
   * @ingroup Contributed
   */
  namespace contrib
  {
    /**
     * @ingroup Contributed
     */
    //! TODO: rewrite to get rid of side-effects and use real RAII
    template<
      class Encapsulation
    >
    class FFT
    {
      public:
        using encap_t = Encapsulation;

      protected:
        struct workspace {
          //fftw_plan           ffft;
          //fftw_plan           ifft;
          //fftw_complex*       wk;
          //complex<typename encap_t::traits::spatial_t>* z;
        };

        map<array<int, Encapsulation::traits::DIM>, shared_ptr<workspace>> workspaces;

      public:
        FFT() = default;
        FFT(const FFT<Encapsulation>& other) = default;
        FFT(FFT<Encapsulation>&& other) = default;
        virtual ~FFT();
        FFT<Encapsulation>& operator=(const FFT<Encapsulation>& other) = default;
        FFT<Encapsulation>& operator=(FFT<Encapsulation>&& other) = default;

        complex<typename Encapsulation::traits::spatial_t>* forward(const shared_ptr<Encapsulation> x);
        void backward(shared_ptr<Encapsulation> x);

        shared_ptr<workspace> get_workspace(const array<int, Encapsulation::traits::DIM>& ndofs);
    };
  }  // ::pfasst::contrib
}  // ::pfasst

#include "pfasst/contrib/fft_impl.hpp"

#endif  // _PFASST__CONTRIB__FFT_HPP_
