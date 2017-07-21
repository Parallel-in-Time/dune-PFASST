// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PFASST_FUNCTIONAL_ENERGYFUNCTIONALCONSTRAINEDLINEARIZATION_HH
#define DUNE_PFASST_FUNCTIONAL_ENERGYFUNCTIONALCONSTRAINEDLINEARIZATION_HH


#include <cstddef>
#include <dune/common/hybridutilities.hh>



namespace Dune {
namespace TNNMG {



/**
 * \brief A constrained linearization for EnergyFunctional
 *
 */
template<class F, class BV>
class EnergyFunctionalConstrainedLinearization
{
public:
    using Matrix = typename F::Matrix;
    using Vector = typename F::Vector;
    using BitVector = BV;
    using ConstrainedVector = Vector;

    using This = EnergyFunctionalConstrainedLinearization<F, BV>;

    EnergyFunctionalConstrainedLinearization(const F& f, const BitVector& ignore) :
        f_(f),
        ignore_(ignore),
        truncationTolerance_(1e-10)
    {}

    void bind(const Vector& x)
    {
        negativeGradient_ = derivative(f_)(x);
        negativeGradient_ *= -1;
        hessian_= hessianMatrix(f_)(x);
        truncationFlags_ = ignore_;
    }

    template<class NV, class NBV>
    static void truncateVector(NV& x, const NBV& truncationFlags)
    {
      namespace H = Dune::Hybrid;
      H::ifElse(IsNumber<NV>(), [&](auto id){
          if (id(truncationFlags))
          id(x) = 0;
          }, [&](auto id){
          H::forEach(H::integralRange(H::size(id(x))), [&](auto&& i) {
              This::truncateVector(x[i], truncationFlags[i]);
              });
          });
    }


    void extendCorrection(ConstrainedVector& cv, Vector& v) const
    {
      v=cv;
      truncateVector(v, truncationFlags_);
    }

    const BitVector& truncated() const
    {
        return truncationFlags_;
    }

    const auto& negativeGradient() const
    {
        return negativeGradient_;
    }

    const auto& hessian() const
    {
        return hessian_;
    }

private:
    const F& f_;
    const BitVector& ignore_;

    double truncationTolerance_;

    Vector negativeGradient_;
    Matrix hessian_;
    BitVector truncationFlags_;
};



} // end namespace TNNMG
} // end namespace Dune




#endif // DUNE_TNNMG_FUNCTIONAL_ENERGYFUNCTIONALCONSTRAINEDLINEARIZATION
