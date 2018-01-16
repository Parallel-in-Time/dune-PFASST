#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

// include the operator from tutorial 01
#include"nonlinearpoissonfem.hh"

/** a local operator for solving the equation
 *
 *  \partial_t u - \Delta u + q(u) = f   in \Omega
 *                               u = g   on \Gamma_D\subseteq\partial\Omega
 *              - \nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *                               u = u_0 at t=t_0
 *
 * (spatial part!) with conforming finite elements on all types of grids in any dimension
 *
 * \tparam Param    A parameter class providing all coefficient functions in the PDE
 * \tparam FEM      Type of a finite element map
 */
template<typename Param, typename FEM>
class NonlinearHeatFEM :
  public NonlinearPoissonFEM<Param,FEM>,
  public Dune::PDELab::
    InstationaryLocalOperatorDefaultMethods<double>
{
  Param& param;
public:
  //! Pass on constructor
  NonlinearHeatFEM (Param& param_, int incrementorder_=0)
    : NonlinearPoissonFEM<Param,FEM>(param_,incrementorder_),
    param(param_)
  {}
  //! set time for subsequent evaluation
  void setTime (double t) {
    param.setTime(t);
  }
};

/** a local operator for the mass operator (L_2 integral)
 *
 * \f{align*}{
 \int_\Omega uv dx
 * \f}
 * \tparam FEM      Type of a finite element map
 */
template<typename FEM>
class L2
  : public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::
  InstationaryLocalOperatorDefaultMethods<double>
{
  typedef typename FEM::Traits::FiniteElementType::
     Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    // types & dimension
    typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsu.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions
        auto& phihat =
          cache.evaluateFunction(ip.position(),
                                 lfsu.finiteElement().localBasis());

        // evaluate u
        RF u=0.0;
        for (size_t i=0; i<lfsu.size(); i++) u += x(lfsu,i)*phihat[i];

        // integrate u*phi_i
        RF factor = ip.weight() * geo.integrationElement(ip.position());
        for (size_t i=0; i<lfsv.size(); i++)
          r.accumulate(lfsv,i,u*phihat[i]*factor);
      }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    // types & dimension
    typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsu.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions
        auto& phihat =
          cache.evaluateFunction(ip.position(),
                                 lfsu.finiteElement().localBasis());

        // integrate phi_j*phi_i
        RF factor = ip.weight() * geo.integrationElement(ip.position());
        for (size_t j=0; j<lfsu.size(); j++)
          for (size_t i=0; i<lfsv.size(); i++)
            mat.accumulate(lfsv,i,lfsu,j,phihat[j]*phihat[i]*factor);
      }
  }

  //! apply local jacobian of the volume term -> nonlinear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,x,lfsv,r);
  }

};
