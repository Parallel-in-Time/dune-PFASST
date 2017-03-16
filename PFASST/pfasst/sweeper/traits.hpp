#ifndef _PFASST__SWEEPER__TRAITS_HPP_
#define _PFASST__SWEEPER__TRAITS_HPP_

#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  /**
   * Type Traits for the Sweepers.
   *
   * @tparam EncapsulationTraits  Traits specifying the data type for the sweeper using this
   *                              sweeper_traits; for supported types, see encap::encap_traits
   * @tparam Ts                   further optional parameters potentially used by specializations
   *
   * @ingroup Traits
   */
  template<
    class EncapsulationTraits,
    class... Ts
  >
  struct sweeper_traits
  {
    //! type of the Encapsulation traits
    using encap_traits = EncapsulationTraits;
    //! type of the user data encapsulation
    using encap_t = encap::Encapsulation<EncapsulationTraits>;
    //! type of the temporal domain
    using time_t = typename encap_traits::time_t;
    //! type of the spacial domain
    using spatial_t = typename encap_traits::spatial_t;
  };
}  // ::pfasst

#endif  // _PFASST__SWEEPER__TRAITS_HPP_
