#ifndef _PFASST__SWEEPER__IMEX_HPP_
#define _PFASST__SWEEPER__IMEX_HPP_

#include <memory>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/sweeper/sweeper.hpp"

//#include "../../../src/finite_element_stuff/fe_manager_fp.hpp"

namespace pfasst
{
  /**
   * Implicit-Explicit Spectral Deferred Corrections Sweeper.
   *
   * The IMEX SDC sweeper uses a specific splitting of the problem equation's right hand side
   * according to the way the parts are solved: either explicit or implicit.
   *
   * Assuming the problem equation @f$ \frac{\partial \vec{u}}{\partial t} = F(\vec{u}, t) @f$,
   * the splitting is defined by @f$ F(\vec{u},t) = F_E(\vec{u},t) + F_I(\vec{u},t) @f$.
   *
   * @f[
   *   \vec{u}_{n+1} = \vec{u}_n + \Delta_{t_{n+1}} \left( F_E(\vec{u}_{n+1},t_{n+1})
   *    - F_E(\vec{u}_n, t_n) + F_I(\vec{u}_{n+1}, t_{n+1}) - F_I(\vec{u}_n,t_n) + Q F_n \right)
   * @f]
   *
   * @ingroup Sweepers
   */
  template<
    class SweeperTrait,
    class BaseFunction,
    typename Enabled = void
  >
  class IMEX
    : public Sweeper<SweeperTrait, BaseFunction, Enabled>
  {
    public:
      //! @copydoc Sweeper::traits
      using traits = SweeperTrait;
      	
    //protected:


      //typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType;


      //typedef Dune::BlockVector<Dune::FieldVector<double,1> > VectorType;



      //std::shared_ptr<fe_manager> FinEl;

      //std::shared_ptr<MatrixType> M_dune;

      bool is_coarse;
      
      bool dontmatter= true;
      
      //typename SweeperTrait::encap_traits::mass_t M_Type;	
      //std::shared_ptr<MatrixType> M_dune;
      //std::shared_ptr<MatrixType> A_dune;

      //typename SweeperTrait::encap_traits::mass_t M_dune;
      //typename SweeperTrait::encap_traits::mass_t A_dune;
      //MatrixType A_dune;
      typename SweeperTrait::encap_traits::mass_t M_dune;
      typename SweeperTrait::encap_traits::mass_t A_dune;	
      
      //Dune::BCRSMatrix <Dune::FieldMatrix<double, 2, 2>> M_dune;

      /**
       * Delta Matrix for the explicit function evaluations.
       *
       * This is a lower triangonal matrix.
       *
       * @see IMEX::compute_delta_matrices()
       */
      Matrix<typename traits::time_t> _q_delta_expl;
      /**
       * Delta Matrix for the implicit function evaluations.
       *
       * This is a lower triangonal matrix with non-zero diagonal.
       *
       * @see IMEX::compute_delta_matrices()
       */
      Matrix<typename traits::time_t> _q_delta_impl;

      //! Cache for the integral @f$ QF_n @f$.
      vector<shared_ptr<typename traits::encap_t>> _q_integrals;
      //! Cache for the explicit function evaluations @f$ F_E(\vec{u}_n, t_n) @f$.
      //vector<shared_ptr<typename traits::encap_t>> _expl_rhs;
      //! Cache for the implicit function evaluations @f$ F_I(\vec{u}_n, t_n) @f$.
      vector<shared_ptr<typename traits::encap_t>> _impl_rhs;
      vector<shared_ptr<typename traits::encap_t>> _impl_rhs_restrict;
      

      //! Counter for total number of @f$ F_E(\vec{u},t) @f$ evaluations.
      //size_t _num_expl_f_evals;
      //! Counter for total number of @f$ F_I(\vec{u},t) @f$ evaluations.
      size_t _num_impl_f_evals;
      //! Counter for total number of implicit solves.
      size_t _num_impl_solves;

      //! @{
      /**
       * @copybrief Sweeper::integrate_end_state
       */
      virtual void integrate_end_state(const typename SweeperTrait::time_t& dt) override;
      /**
       * @copybrief Sweeper::compute_residuals(const bool&)
       */
      virtual void compute_residuals(const bool& only_last) override;

      //virtual void test_m(MatrixType M) override;
      /**
       * @copybrief Sweeper::compute_residuals()
       */
      virtual void compute_residuals() override;
      /**
       * @copybrief Sweeper::initialize()
       */
      virtual void initialize() override;
      //! @}

      //! @name Problem Equation Evaluation
      //! @{
      /**
       * Evaluation of the explicit part of the problem equation's right hand side: @f$ F_E(\vec{u}, t) @f$.
       *
       * @param[in] t  time point to evaluate @f$ F_E(\vec{u},t) @f$ at
       * @param[in] u  spatial data to be passed to @f$ F_E(\vec{u},t) @f$
       * @returns spatial values of function evaluation
       */
      //virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
      //                                                                     const shared_ptr<typename SweeperTrait::encap_t> u);
      
      /**
       * Evaluation of the implicit part of the problem equation's right hand side: @f$ F_I(\vec{u}, t) @f$.
       *
       * @param[in] t  time point to evaluate @f$ F_I(\vec{u},t) @f$ at
       * @param[in] u  spatial data to be passed to @f$ F_I(\vec{u},t) @f$
       * @returns spatial values of function evaluation
       */
      /*virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                                           const shared_ptr<typename SweeperTrait::encap_t> u);*/


      virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                                           const shared_ptr<typename SweeperTrait::encap_t> u);


      /**
       * Implicitly solving the implicit SDC equation.
       *
       * The implicit SDC equation
       * @f[
       *   \vec{u}_{n+1} - \Delta_{t_{n+1}} F_I(\vec{u}_{n+1}, t_{n+1}) = \vec{u}_n
       *    + \Delta_{t_{n+1}} \left( F_E(\vec{u}_{n+1},t_{n+1}) - F_E(\vec{u}_n, t_n)
       *    - F_I(\vec{u}_n,t_n) + Q F_n \right)
       * @f]
       * is solved for @f$ \vec{u}_{n+1} @f$ and @f$ F_I(\vec{u}_{n+1}, t_{n+1}) @f$.
       *
       * @param[in,out] f    values of implicit function evaluation @f$ F_I(\vec{u}_{n+1}, t_{n+1}) @f$
       * @param[in,out] u    desired spatial data @f$ \vec{u}_{n+1} @f$
       * @param[in]     t    _target_ time point @f$ t_{n+1} @f$
       * @param[in]     dt   width of the current time interval
       * @param[in]     rhs  right hand side w.r.t. to SDC as defined by equation above
       */
      virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                  shared_ptr<typename SweeperTrait::encap_t> u,
                                  const typename SweeperTrait::time_t& t,
                                  const typename SweeperTrait::time_t& dt,
                                  const shared_ptr<typename SweeperTrait::encap_t> rhs);
      //! @}

      //! @{
      /**
       * Computes @f$ Q_\Delta @f$ matrices for IMEX scheme.
       */
      virtual void compute_delta_matrices();
      //! @}

    public:
      bool                                            _write = false;
  
      //! @{
      explicit IMEX();
      IMEX(const IMEX<SweeperTrait, BaseFunction, Enabled>& other) = default;
      IMEX(IMEX<SweeperTrait, BaseFunction, Enabled>&& other) = default;
      virtual ~IMEX() = default;
      IMEX<SweeperTrait, BaseFunction, Enabled>& operator=(const IMEX<SweeperTrait, BaseFunction, Enabled>& other) = default;
      IMEX<SweeperTrait, BaseFunction, Enabled>& operator=(IMEX<SweeperTrait, BaseFunction, Enabled>&& other) = default;
      //! @}

      //! @name Configuration and Setup
      //! @{
      /**
       * @copybrief Sweeper::setup()
       *
       * Doesn't do anything special beside calling `Sweeper::setup()`.
       */
      virtual void setup() override;
      //! @}

      //! @name Prediction Step
      //! @{
      /**
       * @copybrief Sweeper::pre_predict()
       *
       * Doesn't do anything special beside calling `Sweeper::pre_predict()`.
       */
      virtual void pre_predict() override;
      /**
       * @copybrief Sweeper::predict()
       *
       * Uses a first order Euler method, which is coincidentally the same than a SDC sweep with
       * the integral @f$ Q F_n @f$ set to zero:
       * @f[
       *   \vec{u}_{n+1} = \vec{u}_n + \Delta_{t_{n+1}} \left( F_E(\vec{u}_{n+1},t_{n+1})
       *    - F_E(\vec{u}_n, t_n) + F_I(\vec{u}_{n+1}, t_{n+1}) - F_I(\vec{u}_n,t_n) \right)
       * @f]
       */
      virtual void predict() override;
      /**
       * @copybrief Sweeper::post_predict()
       *
       * Doesn't do anything special beside calling `Sweeper::post_predict()`.
       */
      virtual void post_predict() override;
      //! @}

      //! @name Sweep Step
      //! @{
      /**
       * @copybrief Sweeper::pre_sweep()
       *
       * Computes the integral @f$ Q F_n @f$ for this iteration.
       */
      virtual void pre_sweep() override;
      /**
       * @copybrief Sweeper::sweep()
       *
       * Propagates the solution node-wise via the SDC equation.
       */
      virtual void sweep() override;
      /**
       * @copybrief Sweeper::post_sweep()
       *
       * Doesn't do anything special beside calling `Sweeper::post_sweep()`.
       */
      virtual void post_sweep() override;
      //! @}

      //! @name Post-Step and Utilities
      //! @{
      /**
       * @copybrief Sweeper::post_step()
       *
       * Doesn't do anything special beside calling `Sweeper::post_step()`.
       */
      virtual void post_step() override;
      /**
       * @copybrief Sweeper::advance(const size_t&)
       * @copydoc Sweeper::advance(const size_t&)
       *
       * Replaces the initial state with the current end state.
       */
      virtual void advance(const size_t& num_steps) override;
      /**
       * @copybrief Sweeper::advance()
       * @copydoc Sweeper::advance()
       */
      virtual void advance() override;
      /**
       * @copybrief Sweeper::reevaluate(const bool)
       *
       * For each node calls `evaluate_rhs_expl()` and `evaluate_rhs_impl()`.
       *
       * @param[in] initial_only  if `true` only evaluates at the initial state
       */
      virtual void reevaluate(const bool initial_only) override;
      /**
       * @copybrief Sweeper::reevaluate()
       * @copydoc Sweeper::reevaluate()
       */
      virtual void reevaluate() override;
      /**
       * @copybrief Sweeper::integrate
       *
       * @param[in] dt  width of the time interval
       */
      virtual vector<shared_ptr<typename SweeperTrait::encap_t>> integrate(const typename SweeperTrait::time_t& dt) override;
      //! @}
            
      virtual vector<shared_ptr<typename SweeperTrait::encap_t>> integrate_new(const typename SweeperTrait::time_t& dt);
  
      //shared_ptr<typename SweeperTrait::encap_traits::mass_t> get_M_dune(){return this->M_dune;};
      shared_ptr<MatrixType> get_M_dune() {return make_shared<MatrixType>(M_dune);};
  };
}  // ::pfasst

#include "pfasst/sweeper/FE_impl_impl.hpp"

#endif  // _PFASST__SWEEPER__IMEX_HPP_
