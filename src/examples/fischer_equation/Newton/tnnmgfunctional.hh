// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=8 sw=2 et sts=2:
#ifndef DUNE_PFASST_FUNCTIONALS_ENERGYFUNCTIONAL_HH
#define DUNE_PFASST_FUNCTIONALS_ENERGYFUNCTIONAL_HH



#include <dune/common/concept.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/solvers/common/resize.hh>
#include <dune/solvers/common/algorithm.hh>
#include <dune/solvers/common/copyorreference.hh>
#include <dune/solvers/common/interval.hh>

#include <dune/tnnmg/iterationsteps/nonlineargsstep.hh>


namespace Dune {
namespace TNNMG {

namespace ComponentMaps {
  /* Maps that act on the numeric components of (possibly non-homogenous or static) containers.
   */

  /** Apply a function f componentwise to a vector v and store it in u, i.e.
   * u_i = f(v_i)
   */
  template<class F, class V, class V2>
  void apply_componentwise(V&& u, V2&& v, F&& f) {
    namespace H = Dune::Hybrid;
    H::ifElse(Dune::IsNumber<std::decay_t<V>>(), [&](auto&& id) {
        v= f(id(u));},
        [&](auto&& id) {
          H::forEach(H::integralRange(H::size(id(v))), [&](auto i) {
            apply_componentwise(H::elementAt(v,i), H::elementAt(u,i), f);
            });
        });
  }

  /** Multiply the i-th component of u by the i-th component of v, i.e.
   * u_i <- u_i*v_i
   */
  template<class V, class V2>
  void multiply_componentwise(V&& u, V2&& v) {
    namespace H = Dune::Hybrid;
    H::ifElse(Dune::IsNumber<std::decay_t<V>>(), [&](auto&& id) {
        v*=id(u);},
        [&](auto&& id) {
          H::forEach(H::integralRange(H::size(id(v))), [&](auto i) {
            multiply_componentwise(H::elementAt(v,i), H::elementAt(u,i));
            });
        });
  }

  /** Write a vector to the diagonal of a matrix 
   */
  template<class M, class V>
  void write_diagonal(M& m, V&& v) {
    namespace H = Dune::Hybrid;
    H::ifElse(Dune::IsNumber<std::decay_t<V>>(), [&](auto&& id) {
        m=id(v);},
        [&](auto&& id) {
          H::forEach(H::integralRange(H::size(id(v))), [&](auto i) {
            write_diagonal(H::elementAt(H::elementAt(m,i), i), H::elementAt(v,i));
            });
        });
  }

  /** Sum all entries of v to sum */
  template<class V, class Field>
  void sum_entries(V&& v, Field* sum) {
    namespace H = Dune::Hybrid;
    H::ifElse(Dune::IsNumber<std::decay_t<V>>(), [&](auto&& id) {
        (*sum)+=id(v);},
        [&](auto&& id) {
          H::forEach(H::integralRange(H::size(id(v))), [&](auto i) {
            sum_entries(H::elementAt(v,i), sum);
            });
        });
  }
}


template<class M, class V, class PHI, class PHIPRIME, class PHI2PRIME, class R=double>
class EnergyFunctional;

template<class M, class V, class PHIPRIME, class R=double>
class EnergyDirectionalRestriction;

/** \brief Shifted Energy functional
 *
 *  \tparam M Global matrix type
 *  \tparam V global vector type
 *  \tparam R Range type
 *  \tparam PHI Type of phi
 *  \tparam PHIPRIME Type of phi'
 *  \tparam PHI2PRIME Type of phi''
 */
template<class M, class V, class PHI, class PHIPRIME, class PHI2PRIME, class R>
class ShiftedEnergyFunctional
{
public:

  using Matrix = M;
  using Vector = V;

  using Range = R;

  ShiftedEnergyFunctional(const Matrix& quadraticPart, const Vector& linearPart, const Vector& weights, const Vector& origin, const Vector& originalOrigin, const PHI& phi,
      const PHIPRIME& dphi, const PHI2PRIME& d2phi) :
    quadraticPart_(&quadraticPart),
    originalLinearPart_(&linearPart),
    origin_(&origin),
    originalOrigin_(&originalOrigin),
    weights_(&weights),
    phi_(phi),
    dphi_(dphi),
    d2phi_(d2phi)
  {}

  // J(v)
  Range operator()(const Vector& v) const
  {
    DUNE_THROW(Dune::NotImplemented, "no evaluation yet!");
  }

  Vector& origin()
  {
    return *origin_;
  }

  const Vector& origin() const
  {
    return *origin_;
  }

  const Vector& weights() const {
    return *weights_;
  }

  void updateOrigin()
  {}

  void updateOrigin(std::size_t i)
  {}

  const Matrix& quadraticPart() const
  {
    return *quadraticPart_;
  }

  const Vector& originalLinearPart() const
  {
    return *originalLinearPart_;
  }

  const Vector& originalOrigin() const
  {
    return *originalOrigin_;
  }

  const auto getPhi() const {
    return phi_;
  }

  const auto getPhiPrime() const {
    return dphi_;
  }

  const auto getPhi2Prime() const {
    return d2phi_;
  }

protected:
  const Matrix* quadraticPart_;
  const Vector* originalLinearPart_;
  const Vector* origin_;
  const Vector* originalOrigin_;
  const Vector* weights_;
  PHI phi_;
  PHIPRIME dphi_;
  PHI2PRIME d2phi_;
};

/** Restrict a ShiftedEnergyFunctional to the coordinate direction of the i-th block */
template<class M, class V, class PHI, class PHIPRIME, class PHI2PRIME, class R, class Index>
auto coordinateRestriction(const ShiftedEnergyFunctional<M,V,PHI, PHIPRIME, PHI2PRIME,R>& f, const Index& i) {
  using LocalMatrix = std::decay_t<decltype(f.quadraticPart()[i][i])>;
  using LocalVector = std::decay_t<decltype(f.originalLinearPart()[i])>;

  using namespace Dune::Solvers;
  namespace H = Dune::Hybrid;

  LocalVector ri = f.originalLinearPart()[i];
  const LocalMatrix* Aii_p = nullptr;

  const auto& Ai = f.quadraticPart()[i];
  sparseRangeFor(Ai, [&](auto&& Aij, auto&& j) {
    // TODO Here we must implement a wrapper to guarantee that this will work with proxy matrices!
    H::ifElse(H::equals(j, i), [&](auto&& id){
      Aii_p = id(&Aij);});
    Imp::mmv(Aij, f.origin()[j], ri, PriorityTag<1>()); // ri -= Aij*f..[j]
  });


  return EnergyFunctional<LocalMatrix, LocalVector, PHI, PHIPRIME, PHI2PRIME, R>(*Aii_p, std::move(ri), f.weights()[i], f.getPhi(), f.getPhiPrime(), f.getPhi2Prime(), f.originalOrigin()[i] +f.origin()[i]);
}



/** \brief Directional restriction of the EnergyFunctional
 *
 *  \tparam M Global matrix type
 *  \tparam V global vector type
 *  \tparam PHIPRIME Type of phi'
 *  \tparam R Range type
 */
template<class M, class V, class PHIPRIME, class R>
class EnergyDirectionalRestriction
{
public:

  using GlobalMatrix = M;
  using GlobalVector = V;

  using Matrix = typename GlobalVector::field_type;
  using Vector = typename GlobalVector::field_type;

  using Range = R;

  EnergyDirectionalRestriction(const GlobalMatrix& matrix, const GlobalVector& linearTerm, const GlobalVector& weights, const GlobalVector& origin, const GlobalVector& direction, PHIPRIME dphi) :
    origin_(origin),
    direction_(direction),
    weights_(weights),
    dphi_(dphi)
  {
    GlobalVector temp = linearTerm;
    matrix.mv(direction, temp); // temp=Ax

    quadraticPart_ = temp*direction; //quadr: <Ax,x>
    linearPart_ = -(linearTerm*direction) + temp*origin; //-<b,x>+ <Ax,u>
  }
  EnergyDirectionalRestriction() {
  }

  Range operator()(const Vector& v) const
  {
    DUNE_THROW(Dune::NotImplemented, "Evaluation of EnergyDirectionalRestriction not implemented");
  }

  const Matrix& quadraticPart() const
  {
    return quadraticPart_;
  }

  const Vector& linearPart() const
  {
    return linearPart_;
  }

  template<class T>
  Dune::Solvers::Interval<Range> subDifferential(T&& t) const {

    // compute derivative at t. The "sub"differential is a scalar as our function is diff'able
    auto derivative = linearPart_ + t*quadraticPart_;

    // add nonlinear part
    auto tmp1 = direction_.get();
    tmp1*=t;
    tmp1+=origin_.get();
    auto tmp2= tmp1;
    ComponentMaps::apply_componentwise(tmp1, tmp2, dphi_); // tmp2 = phi(tmp1) (component-wise)
    // .* x and .*w
    ComponentMaps::multiply_componentwise(direction_.get(), tmp2);
    ComponentMaps::multiply_componentwise(weights_.get(), tmp2);
    // sum all entries on top of derivative
    ComponentMaps::sum_entries(tmp2, &derivative);

    auto I = Dune::Solvers::Interval<Range>(derivative);

    return I;
  }

  // As our functional is not constrained, the domain is the whole real line
  Dune::Solvers::Interval<Range> domain() const {
    return Dune::Solvers::Interval<Range>{std::numeric_limits<Range>::lowest(), std::numeric_limits<Range>::max()};
  }

protected:

  Solvers::CopyOrReference<GlobalVector> origin_;
  Solvers::CopyOrReference<GlobalVector> direction_;
  Solvers::CopyOrReference<GlobalVector> weights_;
  PHIPRIME dphi_;
  Matrix quadraticPart_;
  Vector linearPart_;
};



/** \brief An EnergyFunctional with smooth potentials that act componentwise
 *
 *  \tparam M Matrix type
 *  \tparam V Vector type
 *  \tparam PHI Type of the scalar map phi
 *  \tparam PHIPRIME Type of the scalar map phi'
 *  \tparam PHI2PRIME Type of the scalar map phi''
 *  \tparam R Range type
 */
template<class M, class V,class PHI, class PHIPRIME, class PHI2PRIME, class R>
class EnergyFunctional
{
public:

  using Matrix = std::decay_t<M>;
  using Vector = std::decay_t<V>;
  using Range = R;

  // constructor:
  EnergyFunctional(const M& matrix, const V& b, const V& w, const PHI& phi, const PHIPRIME& phiprime, const PHI2PRIME& phi2prime, const V& offset) :
    matrix_(matrix),
    b_(b),
    w_(w),
    offset_(offset),
    phi_(phi),
    dphi_(phiprime),
    d2phi_(phi2prime)
  {}

  EnergyFunctional(const M& matrix, const V& b, const V& w, const PHI& phi, const PHIPRIME& phiprime, const PHI2PRIME& phi2prime) :
    EnergyFunctional(matrix,b,w,phi,phiprime,phi2prime, V(b.size())) {}

  // F(v)
  Range operator()(const Vector& v) const
  {
    DUNE_THROW(Dune::NotImplemented, "no evaluation yet");
  }

  /**
   * \brief Access derivative
   */
  friend auto derivative(const EnergyFunctional& func)
  {
    return [&](auto&& uu) {// Au -b +phi'(u).*w
      auto tmp = uu;
      func.matrix_.get().mv(uu,tmp); // tmp = Au
      tmp-=func.b_.get();
      auto tmp2 = uu;

      ComponentMaps::apply_componentwise(uu, tmp2, func.dphi_);
      ComponentMaps::multiply_componentwise(func.w_.get(), tmp2);

      tmp+=tmp2;
      return tmp;
    };
  }

  /**
   * \brief Compute the Hessian
   *
   * i.e. compute A + diag(phi''(u_i)*w_i)
   */
  friend auto hessianMatrix(const EnergyFunctional& func) {
    return [&](auto&& uu) {
      auto tmp = func.matrix_.get();
      using Field = typename decltype(tmp)::field_type;
      auto tmp_matrix = func.matrix_.get();
      tmp_matrix=0;

      auto tmp_vector = uu;
      tmp_vector = 0;

      ComponentMaps::apply_componentwise(uu, tmp_vector, func.d2phi_); // tmp_vector = d2phi(u)
      ComponentMaps::multiply_componentwise(func.w_.get(), tmp_vector); // tmp_vector .* w_
      ComponentMaps::write_diagonal(tmp_matrix, tmp_vector); // write entries on the matrix' diagonal

      tmp+=tmp_matrix;
      return tmp;
    };
  }

  friend auto directionalRestriction(const EnergyFunctional& f, const Vector& origin, const Vector& direction)
    -> EnergyDirectionalRestriction<Matrix, Vector, PHIPRIME, Range>
  {
    return EnergyDirectionalRestriction<Matrix, Vector, PHIPRIME,Range>(f.matrix_.get(), f.b_.get(), f.w_.get(), origin, direction, f.dphi_);
  }

  friend auto shift(const EnergyFunctional& f, const Vector& origin)
  {
    return ShiftedEnergyFunctional<Matrix, Vector, PHI, PHIPRIME, PHI2PRIME, Range>(f.matrix_.get(), f.b_.get(), f.w_.get(), origin, f.offset_, f.phi_, f.dphi_, f.d2phi_);
  }

  // As our functional is not constrained, the domain is the whole real line
  Dune::Solvers::Interval<Range> domain() const {
    return Dune::Solvers::Interval<Range>{std::numeric_limits<Range>::lowest(), std::numeric_limits<Range>::max()};
  }

  template<class T>
  Dune::Solvers::Interval<Range> subDifferential(T&& t) const {

    // compute derivative at t. The "sub"differential is a scalar as our function is diff'able
    auto sum = 0.0;
    namespace H = Dune::Hybrid;
    H::ifElse(Dune::IsNumber<std::decay_t<Vector>>::value, [&](auto&& id) {
      sum = -b_.get() + t*matrix_.get();
      sum += dphi_(offset_ + t)*w_.get();
      }
      ,[](auto) {DUNE_THROW(Dune::Exception, "This case should not happen, your nonlinear GS is probably not nested");}
    );
    auto I = Dune::Solvers::Interval<Range>(sum);

    return I;
  }

protected:

  Solvers::ConstCopyOrReference<M> matrix_;
  Solvers::ConstCopyOrReference<V> w_;
  Solvers::ConstCopyOrReference<V> b_;
  V offset_;
  PHI phi_;
  PHIPRIME dphi_;
  PHI2PRIME d2phi_;
};




} // end namespace PFASST
} // end namespace Dune

#endif // DUNE_PFASST_FUNCTIONALS_ENERGYFUNCTIONAL_HH
