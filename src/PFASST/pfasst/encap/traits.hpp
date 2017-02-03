#ifndef _PFASST__ENCAP__TRAITS_HPP_
#define _PFASST__ENCAP__TRAITS_HPP_

#include <type_traits>
#include <vector>
using std::vector;

#include "pfasst/globals.hpp"


namespace pfasst
{
  namespace encap
  {
    /**
     * @ingroup Tags
     */
    struct encap_data_tag
    {};

    /**
     * @ingroup Tags
     */
    struct vector_encap_tag
      : public encap_data_tag
    {};

    /**
     * @ingroup Tags
     */
    struct eigen3_encap_tag
      : public encap_data_tag
    {};

    /**
     * Type Traits for encapsulation of generic user data types.
     *
     * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
     * @tparam SpatialPrecision the spatial data precision
     * @tparam Dim              number of spatial dimensions
     * @tparam DataT            the actual data type encapsulated
     * @tparam Ts               further optional parameters potentially used by specializations
     *
     * @ingroup Traits
     */
    template<
      class TimePrecision,
      class SpatialPrecision,
      size_t Dim,
      class DataT,
      class... Ts
    >
    struct encap_traits
    {
      //! type of the time precision, e.g. precision of integration nodes
      using time_t = TimePrecision;

      //! type of the spatial precision
      using spatial_t = SpatialPrecision;

      //! type of the encapsulated data
      using data_t = DataT;

      using tag_t = encap_data_tag;

      using dim_t = std::integral_constant<size_t, Dim>;

      //! number of spatial dimensions
      static constexpr size_t DIM = Dim;
    };


    /**
     * Specialized Type Traits for encapsulation of std::vector.
     *
     * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
     * @tparam SpatialPrecision the spatial data precision
     *
     * @ingroup Traits
     */
    template<
      class TimePrecision,
      class SpatialPrecision,
      size_t Dim
    >
    struct vector_encap_traits
      : public encap_traits<TimePrecision, SpatialPrecision, Dim, vector<SpatialPrecision>>
    {
      using time_t = TimePrecision;
      using spatial_t = SpatialPrecision;
      using data_t = vector<spatial_t>;
      using tag_t = vector_encap_tag;
      using dim_t = std::integral_constant<size_t, Dim>;
      static constexpr size_t DIM = Dim;
    };


    /**
     * Specialized Type Traits for encapsulation of std::vector.
     *
     * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
     * @tparam SpatialPrecision the spatial data precision
     *
     * @ingroup Traits
     */
    template<
      class TimePrecision,
      class SpatialPrecision,
      size_t Dim
    >
    struct eigen3_encap_traits
      : public encap_traits<TimePrecision, SpatialPrecision, Dim, EigenVector<SpatialPrecision>>
    {
      using time_t = TimePrecision;
      using spatial_t = SpatialPrecision;
      using data_t = EigenVector<spatial_t>;
      using tag_t = eigen3_encap_tag;
      using dim_t = std::integral_constant<size_t, Dim>;
      static constexpr size_t DIM = Dim;
    };
  }  // ::pfasst::encap
}  // ::pfasst

#endif  // _PFASST__ENCAP__TRAITS_HPP_
