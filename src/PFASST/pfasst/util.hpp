#ifndef _PFASST__UTIL_HPP_
#define _PFASST__UTIL_HPP_

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  using std::string;
  using std::vector;

  /**
   * Compares absolute difference against zero.
   *
   * Takes the absolute value of the difference between @p a and @p b and compares it to the machine
   * precision scaled to the magnitude of @p a and @p b up to the desired precision in units in the
   * last place.
   *
   * @tparam    precision  value type of @p a and @p b
   * @param[in] a
   * @param[in] b
   * @param[in] digits     number of digits to compare
   * @ingroup Utilities
   */
  template<typename precision>
  static bool almost_equal(const precision& a, const precision& b,
                           const int digits = std::numeric_limits<precision>::digits)
  {
    using std::abs;

    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return abs(a - b) < std::numeric_limits<precision>::epsilon() * abs(a + b) * digits
    // unless the result is subnormal
           || abs(a - b) < std::numeric_limits<precision>::min();
  }

  /**
   * Compares against smallest machine epsilon.
   *
   * Takes the absolute value of @p a and compares it against `std::numeric_limits::epsilon()`.
   *
   * @tparam    precision  value type of @p a, e.g. `double`
   * @param[in] a          value to compare
   *
   * @returns   `true` if @p a is closer to zero than the machine precision, `false` otherwise.
   *
   * @ingroup Utilities
   */
  template<typename precision>
  static bool almost_zero(const precision& a)
  {
    return std::abs(a) < std::numeric_limits<precision>::epsilon();
  }


  template<
    class... Args,
    typename std::enable_if<sizeof...(Args) == 2>::type* = nullptr
  >
  static constexpr size_t
  linearized_index(const std::tuple<Args...>& indices, const size_t dim_ndofs)
  {
    return std::get<0>(indices) * dim_ndofs + std::get<1>(indices);
  }

  template<
    class... Args,
    typename std::enable_if<sizeof...(Args) == 3>::type* = nullptr
  >
  static constexpr size_t
  linearized_index(const std::tuple<Args...>& indices, const size_t dim_ndofs)
  {
    using std::get;

    return get<0>(indices) * std::pow(dim_ndofs, 2) + get<1>(indices) * dim_ndofs + get<2>(indices);
  }


  template<
    size_t N,
    typename std::enable_if<N == 2>::type* = nullptr
  >
  static const std::tuple<size_t, size_t>
  split_index(const size_t index, const size_t dim_ndofs)
  {
    const size_t xi = index % dim_ndofs;
    const size_t yi = (index - xi) / dim_ndofs;
    return std::make_tuple(yi, xi);
  }

  template<
    size_t N,
    typename std::enable_if<N == 3>::type* = nullptr
  >
  static const std::tuple<size_t, size_t, size_t>
  split_index(const size_t index, const size_t dim_ndofs)
  {
    using std::pow;

    const size_t xi = index % dim_ndofs;
    const size_t zyi = index - xi;
    const size_t zi = (zyi >= pow(dim_ndofs, 2)) ? floor(zyi / pow(dim_ndofs, 2)) : 0;
    const size_t yi = (zyi - zi * pow(dim_ndofs, 2)) / dim_ndofs;
    return std::make_tuple(zi, yi, xi);
  }


  /**
   * A string representation of `std::vector`.
   *
   * Uses `operator<<` on each element of @p vec and seperates each element with @p sep .
   *
   * @tparam    T   type of the vector's elements
   * @param[in] vec vector to print
   * @param[in] sep seperator placed between each element of @p vec
   * @returns string representation of @p vec
   *
   * **Example**:
   * @code{.cpp}
   * std::vector<double> vec{0.0, 0.5, 0.1};
   * join(vec, ", ");  // ==> "[0.0, 0.5, 0.1]"
   * @endcode
   * 
   * @note This funtion is a noop when compiling with `PFASST_NO_LOGGING` defined.
   *
   * @ingroup Utilities
   */
  template<typename T>
  static string join(const vector<T>& vec, const string& sep)
  {
#ifndef PFASST_NO_LOGGING
    std::stringstream os;
    //os << "[" << LOG_FLOAT;
    os << "[" ;
    for (size_t i = 0; i < vec.size() - 1; ++i) {
      os << vec[i] << sep;
    }
    os << vec.back();
    os << "]";
    return os.str();
#else
    UNUSED(vec); UNUSED(sep);
    return "";
#endif
  }
}  // ::pfasst

#endif  // _PFASST__UTIL_HPP_
