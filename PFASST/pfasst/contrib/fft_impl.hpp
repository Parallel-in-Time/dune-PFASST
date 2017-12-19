#include "pfasst/contrib/fft.hpp"

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace contrib
  {
    template<class Encapsulation>
    FFT<Encapsulation>::~FFT()
    {

    }

    template<class Encapsulation>
    complex<typename Encapsulation::traits::spatial_t>*
    FFT<Encapsulation>::forward(const shared_ptr<Encapsulation> x)
    {

    }

    template<class Encapsulation>
    void
    FFT<Encapsulation>::backward(shared_ptr<Encapsulation> x)
    {

    }

    /*template<class Encapsulation>
    shared_ptr<typename FFT<Encapsulation>::workspace>
    FFT<Encapsulation>::get_workspace(const array<int, Encapsulation::traits::DIM>& ndofs)
    {


      return workspaces[ndofs];
    }*/
  }  // ::pfast::contrib
}  // ::pfasst
