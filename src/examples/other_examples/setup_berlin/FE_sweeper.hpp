#ifndef _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_

#include <memory>
#include <type_traits>

using std::shared_ptr;
using std::vector;

#include <vector>

#include <pfasst/sweeper/FE_impl.hpp> //implicit finite element sweeper
#include <pfasst/contrib/fft.hpp>

//#include "fe_manager.hpp"



namespace pfasst
{
  namespace examples
  {
    namespace FE_sweeper
    {
        //dune_sweeper_traits must be somewhere else!
        template<
                class EncapsulationTraits,
                size_t base_order,
                size_t dim,
                class... Ts
        >
        struct dune_sweeper_traits : public sweeper_traits<EncapsulationTraits,  Ts...>
        {
            //! type of the Encapsulation traits
            using encap_traits = EncapsulationTraits;
            //! type of the user data encapsulation
            using encap_t = encap::Encapsulation<EncapsulationTraits>;
            //! type of the temporal domain
            using time_t = typename encap_traits::time_t;
            //! type of the spacial domain
            using spatial_t = typename encap_traits::spatial_t;

            static constexpr size_t BASE_ORDER = base_order;
            static constexpr size_t DIM = dim;
        };

      template<
        class SweeperTrait,
        class BaseFunction,
        typename Enabled = void
      >
      class Heat_FE
        : public IMEX<SweeperTrait, BaseFunction, Enabled>
      {
        static_assert(std::is_same<
                        typename SweeperTrait::encap_t::traits::dim_t,
                        std::integral_constant<size_t, 1>  //ruth_dim
                      >::value,
                      "Heat_FE Sweeper requires 2D data structures");

        public:
          using traits = SweeperTrait;

          static void init_opts();
	  int _iterations{0};
          bool                                            _write = false;

          


        protected:
            

          using spatial_t = typename traits::spatial_t;

 

          size_t nlevel;
	  //std::shared_ptr<fe_manager> FinEl;
          std::shared_ptr<GridType> grid;
	  typedef GridType::LevelGridView GridView;
          std::shared_ptr<BaseFunction> basis;   

          
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
          //explicit Heat_FE(const size_t nelements, const size_t basisorder);
	  //explicit Heat_FE(std::shared_ptr<Dune::Functions::PQkNodalBasis<GridType::LevelGridView,SweeperTrait::BASE_ORDER>> basis, size_t, std::shared_ptr<GridType> grid);
	  explicit Heat_FE(std::shared_ptr<BaseFunction> basis, size_t, std::shared_ptr<GridType> grid);
	  
          
          
          Heat_FE(const Heat_FE<SweeperTrait, BaseFunction, Enabled>& other) = default;
          Heat_FE(Heat_FE<SweeperTrait, BaseFunction, Enabled>&& other) = default;
          virtual ~Heat_FE() = default;
          Heat_FE<SweeperTrait, BaseFunction, Enabled>& operator=(const Heat_FE<SweeperTrait, BaseFunction, Enabled>& other) = default;
          Heat_FE<SweeperTrait, BaseFunction, Enabled>& operator=(Heat_FE<SweeperTrait, BaseFunction, Enabled>&& other) = default;

          virtual void set_options() override;

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);
	  
	  
          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;

          auto get_A_dune() const {
            return this->A_dune;
          }
	  
	
      };
      
      
      
     
      
      
    }  // ::pfasst::examples::heat_FE
  }  // ::pfasst::examples
}  // ::pfasst

#include "FE_sweeper_impl.hpp"





#endif  // _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
