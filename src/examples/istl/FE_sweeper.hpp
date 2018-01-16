#ifndef _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_


#include <memory>
#include <type_traits>

using std::shared_ptr;
using std::vector;

#include <vector>

#include <pfasst/sweeper/FE_impl.hpp>
#include <pfasst/contrib/fft.hpp>

#include "fe_manager.hpp"

//#ifndef PI
//#define PI 3.1415926535897932385
//#endif


namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
        /*template<
                class EncapsulationTraits,
                size_t base_order,
                size_t dim,
                class... Ts
        >
        struct dune_sweeper_traits : public sweeper_traits<EncapsulationTraits,  Ts...>
        {

            using encap_traits = EncapsulationTraits;
            using encap_t = encap::Encapsulation<EncapsulationTraits>;
            using time_t = typename encap_traits::time_t;
            using spatial_t = typename encap_traits::spatial_t;

            //static constexpr size_t BASE_ORDER = base_order;
            //static constexpr size_t DIM = dim;
        };*/

      template<
        class SweeperTrait,
	class Mass,
        typename Enabled = void
      >
      class Heat_FE
        : public IMEX<SweeperTrait, Mass, Enabled>
      {
        /*static_assert(std::is_same<
                        typename SweeperTrait::encap_t::traits::dim_t,
                        std::integral_constant<size_t, 1>
                      >::value,
                      "Heat_FE Sweeper requires 2D data structures");*/

        public:
          using traits = SweeperTrait;

          static void init_opts();

          int finer;

        private:
          using spatial_t = typename traits::spatial_t;

          typename traits::time_t                        _t0{0.0};
          spatial_t                                      _nu{0.1};//0.1
          pfasst::contrib::FFT<typename traits::encap_t> _fft;
          vector<vector<spatial_t>>                      _lap;


          std::shared_ptr<fe_manager> FinEl;
	  


          //MatrixType M_dune;
          //MatrixType A_dune;


          std::shared_ptr<GridType> grid;



          std::shared_ptr<BasisFunction> basis;



        //protected:
      public:  
          /*virtual shared_ptr<typename SweeperTrait::encap_t>
          evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;*/


          shared_ptr<typename SweeperTrait::encap_traits::mass_t> A_dune;	
          virtual shared_ptr<typename SweeperTrait::encap_t>
          evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;

          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                      shared_ptr<typename SweeperTrait::encap_t> u,
                                      const typename SweeperTrait::time_t& t,
                                      const typename SweeperTrait::time_t& dt,
                                      const shared_ptr<typename SweeperTrait::encap_t> rhs) override;

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_error(const typename SweeperTrait::time_t& t);

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                 const typename SweeperTrait::time_t& t);

        public:


          explicit Heat_FE(std::shared_ptr<fe_manager>, size_t);


          Heat_FE(const Heat_FE<SweeperTrait, Mass, Enabled>& other)= default;
          Heat_FE(Heat_FE<SweeperTrait, Mass, Enabled>&& other)= default;
          virtual ~Heat_FE() = default;
          Heat_FE<SweeperTrait, Mass, Enabled>& operator=(const Heat_FE<SweeperTrait, Mass, Enabled>& other) = default;
          Heat_FE<SweeperTrait, Mass, Enabled>& operator=(Heat_FE<SweeperTrait, Mass, Enabled>&& other) = default;

          virtual void set_options() override;

          template<typename Basis>
          void assemble(Basis &basis);

          void assemble();

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);

          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;

          shared_ptr<GridType> get_grid() const;
      };
    }  // ::pfasst::examples::heat_FE
  }  // ::pfasst::examples
}  // ::pfasst

#include "FE_sweeper_impl.hpp"

#endif // _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
