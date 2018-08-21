#ifndef _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP2_
#define _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP2_

#include <memory>
#include <type_traits>

using std::shared_ptr;
using std::vector;

#include <vector>

#include <pfasst/sweeper/FE_impl.hpp>

#include <pfasst/contrib/fft.hpp>
#include "../2d_transfer/fe_manager_fp.hpp"



namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
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
        typename Enabled = void
      >
      class Heat_FE
        : public IMEX<SweeperTrait, Enabled>
      {
        static_assert(std::is_same<
                        typename SweeperTrait::encap_t::traits::dim_t,
                        std::integral_constant<size_t, 1>  //dimension
                      >::value,
                      "Heat_FE Sweeper requires 2D data structures");

        public:
	  int                                             num_solves=0; //how often the linear system is solved 
          double                                         newton = 1e-6;	//tollerance         
          using traits = SweeperTrait;

          static void init_opts();
          int output_level;
          int output=0;
        private:
        

          using spatial_t = typename traits::spatial_t;

          typename traits::time_t                        _t0{0.0};
          double                                     	 _nu{25};
          double                                         _eps{0.4};
          double                                     	 _n{2.0};
          std::shared_ptr<VectorType> w;

	  std::shared_ptr<fe_manager> FinEl;

          std::shared_ptr<GridType> grid;

          std::shared_ptr<BasisFunction> basis;

	  typedef GridType::LevelGridView GridView;

          //using BasisFunction = Dune::Functions::PQkNodalBasis<GridView,1>; //SweeperTrait::BASE_ORDER>;



        protected:


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
	  
	  	  virtual void
	  evaluate_f(shared_ptr<typename SweeperTrait::encap_t> f,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt,
						 const shared_ptr<typename SweeperTrait::encap_t> rhs
						);
						
	  virtual void					
	  evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						);

        public:
          explicit Heat_FE(std::shared_ptr<fe_manager>, size_t);
          //explicit Heat_FE(const size_t nelements, const size_t basisorder);
	      //explicit Heat_FE(std::shared_ptr<Dune::Functions::PQkNodalBasis<GridType::LeafGridView,SweeperTrait::BASE_ORDER>> basis, size_t);
	  
          Heat_FE(const Heat_FE<SweeperTrait, Enabled>& other) = default;
          Heat_FE(Heat_FE<SweeperTrait, Enabled>&& other) = default;
          virtual ~Heat_FE() = default;
          Heat_FE<SweeperTrait, Enabled>& operator=(const Heat_FE<SweeperTrait, Enabled>& other) = default;
          Heat_FE<SweeperTrait, Enabled>& operator=(Heat_FE<SweeperTrait, Enabled>&& other) = default;

          virtual void set_options() override;

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);
	  //virtual shared_ptr<typename SweeperTrait::encap_t> source(const typename SweeperTrait::time_t& t);
	  
	  
          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;

          //shared_ptr<GridType> get_grid() const;
	  
	
      };
    }  // ::pfasst::examples::heat_FE
  }  // ::pfasst::examples
}  // ::pfasst

#include "FE_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
