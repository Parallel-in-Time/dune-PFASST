#ifndef _PFASST__CONTROLLER__STATUS_HPP_
#define _PFASST__CONTROLLER__STATUS_HPP_

#include <limits>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>
using std::shared_ptr;
using std::string;
using std::vector;


#include <better-enums/enum.h>


#ifdef WITH_MPI
  #include <mpi.h>
#endif

#include "pfasst/logging.hpp"


namespace pfasst
{
#ifdef WITH_MPI
  static MPI_Datatype status_data_type;
#endif

  /**
   * @enum PrimaryState::_enumerated
   * @brief Enhanced enumeration for a primary state of the Controller.
   *
   * @note This is an enhanced enumeration using
   *   [Better Enums](http://aantron.github.io/better-enums/index.html).
   *
   * @ingroup Controllers
   */
  ENUM(PrimaryState, int,
    // general state
    CONVERGED         =  0,
    FAILED            =  1,

    // iterating states
    PREDICTING        = 10,
    ITERATING         = 20,

    // intermediate states
    INTER_ITER        = 30,

    UNKNOWN_PRIMARY   = std::numeric_limits<int>::max()
  )

  /**
   * @enum SecondaryState::_enumerated
   * @brief Enhanced enumeration for a secondary state of the Controller.
   *
   * @note This is an enhanced enumeration using
   *   [Better Enums](http://aantron.github.io/better-enums/index.html).
   *
   * @ingroup Controllers
   */
  ENUM(SecondaryState, int,
    // inter-level states
    CYCLE_DOWN        =  1,
    CYCLE_UP          =  2,

    // coarse level
    PRE_ITER_COARSE   =  9,
    ITER_COARSE       = 10,
    POST_ITER_COARSE  = 11,

    // fine level
    PRE_ITER_FINE     = 19,
    ITER_FINE         = 20,
    POST_ITER_FINE    = 21,

    // intermediate states
    CONV_CHECK        = 30,

    UNKNOWN_SECONDARY = std::numeric_limits<int>::max()
  )

  /**
   * @ingroup Utilities
   */
  template<class CharT, class Traits>
  inline std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, const pfasst::PrimaryState& state)
  {
    os << (+state)._to_string();
    return os;
  }

  /**
   * @ingroup Utilities
   */
  template<class CharT, class Traits>
  inline std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, const pfasst::SecondaryState& state)
  {
    os << (+state)._to_string();
    return os;
  }

  // begin of some constexpr stuff
  /**
   * @defgroup StatusDetails Status Details
   * @ingroup Controllers
   * @see The idea for this C++11 constexpr stuff is based on
   *   [this SO answer](http://stackoverflow.com/a/26079954/588243).
   */

  /**
   * Utility typedef for combinations of primary and secondary states.
   *
   * @ingroup StatusDetails
   */
  using PrimarySecondaryMapItem = std::pair<PrimaryState, SecondaryState>;
  /**
   * List of valid combinations of PrimaryState and SecondaryState.
   *
   * @ingroup StatusDetails
   */
  constexpr PrimarySecondaryMapItem VALID_STATES_COMBINATIONS[] = {
    { PrimaryState::CONVERGED,       SecondaryState::UNKNOWN_SECONDARY },

    { PrimaryState::FAILED,          SecondaryState::UNKNOWN_SECONDARY },

    { PrimaryState::PREDICTING,      SecondaryState::UNKNOWN_SECONDARY },
    { PrimaryState::PREDICTING,      SecondaryState::CYCLE_DOWN },
    { PrimaryState::PREDICTING,      SecondaryState::PRE_ITER_COARSE },
    { PrimaryState::PREDICTING,      SecondaryState::ITER_COARSE },
    { PrimaryState::PREDICTING,      SecondaryState::POST_ITER_COARSE },
    { PrimaryState::PREDICTING,      SecondaryState::CYCLE_UP },
    { PrimaryState::PREDICTING,      SecondaryState::PRE_ITER_FINE },
    { PrimaryState::PREDICTING,      SecondaryState::ITER_FINE },
    { PrimaryState::PREDICTING,      SecondaryState::POST_ITER_FINE },

    { PrimaryState::ITERATING,       SecondaryState::UNKNOWN_SECONDARY },
    { PrimaryState::ITERATING,       SecondaryState::CYCLE_DOWN },
    { PrimaryState::ITERATING,       SecondaryState::PRE_ITER_COARSE },
    { PrimaryState::ITERATING,       SecondaryState::ITER_COARSE },
    { PrimaryState::ITERATING,       SecondaryState::POST_ITER_COARSE },
    { PrimaryState::ITERATING,       SecondaryState::CYCLE_UP },
    { PrimaryState::ITERATING,       SecondaryState::PRE_ITER_FINE },
    { PrimaryState::ITERATING,       SecondaryState::ITER_FINE },
    { PrimaryState::ITERATING,       SecondaryState::POST_ITER_FINE },

    { PrimaryState::INTER_ITER,      SecondaryState::UNKNOWN_SECONDARY },
    { PrimaryState::INTER_ITER,      SecondaryState::CONV_CHECK },

    { PrimaryState::UNKNOWN_PRIMARY, SecondaryState::UNKNOWN_SECONDARY }
  };
  /**
   * Total number of valid state combinations as defined in `VALID_STATES_COMBINATIONS`.
   *
   * @ingroup StatusDetails
   */
  constexpr auto VALID_STATES_COMBINATIONS_SIZE = sizeof VALID_STATES_COMBINATIONS / sizeof VALID_STATES_COMBINATIONS[0];

  /**
   * Utility function to check whether two given states represent a valid combination.
   *
   * This is a recursive function iterating over the list of valid combinations, starting with the
   * last one.
   * As soon as both @p primary and @p secondary match, it returns.
   *
   * @param[in] primary    PrimaryState of the combination to check
   * @param[in] secondary  SeconaryState of the combination to check
   * @param[in] range      index in `VALID_STATES_COMBINATIONS` to try
   *
   * @returns `true` when @p primary and @p secondary are a valid combination, `false` otherwise.
   *
   * @ingroup StatusDetails
   */
  static constexpr bool validate_state_combination(PrimaryState primary, SecondaryState secondary,
                                                   int range = VALID_STATES_COMBINATIONS_SIZE) {
      return (range == 0)
               ? false
               : ((VALID_STATES_COMBINATIONS[range - 1].first == primary)
                  ? ((VALID_STATES_COMBINATIONS[range - 1].second == secondary)
                     ? true
                     : validate_state_combination(primary, secondary, range - 1)
                    )
                  : validate_state_combination(primary, secondary, range - 1)
                 );
  }
  //! @}
  // ... end of constexpr stuff


  /**
   * Internal storage for the Status data.
   *
   * @tparam precision numerical type for the time domain, e.g. `double`
   *
   * @note The order of the member variables is actually on purpose to ensure an aligned memory
   *    layout.
   *    This compact memory layout is exploited in the MPI communication routines for the
   *    StatusDetail.
   *
   * @see Status::create_mpi_datatype().
   *
   * @ingroup StatusDetails
   */
  template<
    typename precision
  >
  struct StatusDetail
  {
    //! Primary state.
    PrimaryState   primary_state   = (+pfasst::PrimaryState::UNKNOWN_PRIMARY);
    //! Secondary state.
    SecondaryState secondary_state = (+pfasst::SecondaryState::UNKNOWN_SECONDARY);

    //! Current time step index.
    size_t    step                 = 0;
    //! Total number of time steps to compute.
    size_t    num_steps            = 0;
    //! Current iteration index.
    size_t    iteration            = 0;
    //! Allowed maximum number of iterations per time step.
    size_t    max_iterations       = 0;

    //! Current time step's @f$ t_0 @f$.
    precision time                 = 0.0;
    //! Current time step width.
    precision dt                   = 0.0;
    //! Current time step's @f$ t_{end} @f$.
    precision t_end                = 0.0;
    //! Current maximum of absolute residual norms.
    precision abs_res_norm         = 0.0;
    //! Current maximum of relative residual norms.
    precision rel_res_norm         = 0.0;
  };

  //! @cond NO_DOC
  template<typename precision>
  inline MAKE_LOGGABLE(StatusDetail<precision>, status_detail, os)
  {
    os << "StatusDetail("
       << "t=" << status_detail.time
       << ", dt=" << status_detail.dt
       << ", t_end=" << status_detail.t_end
       << ", step=" << status_detail.step
       << ", num_steps=" << status_detail.num_steps
       << ", iter=" << status_detail.iteration
       << ", iter_max=" << status_detail.max_iterations
       << ", 1st_state=" << (+status_detail.primary_state)._to_string()
       << ", 2nd_state=" << (+status_detail.secondary_state)._to_string()
//       << ", abs_res=" << LOG_FLOAT << status_detail.abs_res_norm
//       << ", rel_res=" << LOG_FLOAT << status_detail.rel_res_norm
       << ", abs_res=" <<  status_detail.abs_res_norm
       << ", rel_res=" <<  status_detail.rel_res_norm
       << ")";
    return os;
  }
  //! @endcond


  /**
   * Status data for Controllers and Sweepers.
   *
   * These Status objects contain all the important runtime information of an executing Controller
   * and Sweeper.
   * They are used inside the Controllers and Sweepers for retrieving information about the current
   * state of the run including the current state, time step and iteration as well as some
   * termination criteria.
   *
   * @tparam precision numerical type of the time domain, e.g. `double`; this is checked at compile
   *                   time to satisfy `std::is_arithmetic`
   *
   * @see StatusDetail for default and initial values
   *
   * @ingroup Controllers
   */
  template<
    typename precision
  >
  class Status
    :   public std::enable_shared_from_this<Status<precision>>
      , public el::Loggable
  {
    static_assert(std::is_arithmetic<precision>::value,
                  "precision type must be arithmetic");

    public:
      //! Numerical type of the time domain.
      using precision_t = precision;

#ifdef WITH_MPI
      //! @name MPI Utilities
      //! @{
      /**
       * Creates a MPI derived data type for StatusDetail.
       *
       * @note You have to call this function exactly once just after calling `pfasst::init()`.
       */
      static inline void create_mpi_datatype();
      /**
       * Destroys/Frees the MPI derived data type for StatusDetail.
       *
       * @note You have to call this function exactly once at the end of your MPI program just
       *   before calling `MPI_Finalize()`.
       */
      static inline void free_mpi_datatype();
      //! @}
#endif

      static_assert(std::is_standard_layout<StatusDetail<precision>>::value,
                    "Status::Detail needs to be have standard layout for MPI derived datatype");

      //! Data storage.
      StatusDetail<precision> _detail;

    public:
      //! @{
      Status() = default;
      Status(const Status<precision>& other) = default;
      Status(Status<precision>&& other) = default;
      virtual ~Status() = default;
      Status<precision>& operator=(const Status<precision>& other) = default;
      Status<precision>& operator=(Status<precision>&& other) = default;
      //! @}

      //! @name Setup
      //! @{
      /**
       * Clears the internal data storage.
       */
      virtual void clear();
      /**
       * Resets the Status for a new time step.
       */
      virtual void reset();
      //! @}

      //! @name Accessors
      //! @{
      /**
       * Accessor for the current time step index.
       */
      virtual size_t& step();
      //! Read-only version of `step()`.
      virtual size_t  get_step() const;

      /**
       * Accessor for the total number of required time steps.
       */
      virtual size_t& num_steps();
      //! Read-only version of `num_steps()`.
      virtual size_t  get_num_steps() const;

      /**
       * Accessor for the current iteration index.
       */
      virtual size_t& iteration();
      //! Read-only version of `iteration()`.
      virtual size_t  get_iteration() const;

      /**
       * Accessor for the total number of allowed iterations.
       */
      virtual size_t& max_iterations();
      //! Read-only version of `max_iterations()`.
      virtual size_t  get_max_iterations() const;

      /**
       * Accessor for the current time step's @f$ t_0 @f$.
       */
      virtual precision& time();
      //! Read-only version of `time()`.
      virtual precision  get_time() const;

      /**
       * Accessor for the current time step's width.
       */
      virtual precision& dt();
      //! Read-only version of `dt()`.
      virtual precision  get_dt() const;

      /**
       * Accessor for the current time step's @f$ t_{end} @f$.
       */
      virtual precision& t_end();
      //! Read-only version of `t_end()`.
      virtual precision  get_t_end() const;

      /**
       * Setter for the primary state.
       *
       * @note This resets the secondary state to `UNKNOWN_SECONDARY` of
       *    `SecondaryState::_enumerated`.
       *
       * @param[in] state new primary state
       */
      virtual void         set_primary_state(const PrimaryState& state);
      //! Read-only accessor for the primary state.
      virtual PrimaryState get_primary_state() const;

      /**
       * Setter for the secondary state.
       *
       * @note This checks the validity of the combination of the current primary state and the new
       *   secondary @p state.
       *
       * @param[in] state new secondary state
       */
      virtual void           set_secondary_state(const SecondaryState& state);
      //! Read-only accessor for the secondary state.
      virtual SecondaryState get_secondary_state() const;

      /**
       * Accessor for the maximum value of the absolute residual norms.
       */
      virtual precision& abs_res_norm();
      //! Read-only version of `abs_res_norm()`.
      virtual precision  get_abs_res_norm() const;

      /**
       * Accessor for the maximum value of the relative residual norms.
       */
      virtual precision& rel_res_norm();
      //! Read-only version of `rel_res_norm()`.
      virtual precision  get_rel_res_norm() const;
      //! @}

      //! @name Communication
      //! @{
      template<class CommT>
      bool probe(shared_ptr<CommT> comm, const int src_rank, const int tag);

      template<class CommT>
      void send(shared_ptr<CommT> comm, const int dest_rank, const int tag, const bool blocking);

      template<class CommT>
      void recv(shared_ptr<CommT> comm, const int src_rank, const int tag, const bool blocking);
      //! @}

      //! @name Logging
      //! @{
      /**
       * Multi-line summary of the Status object.
       *
       * As log files should not include entries spanning multiple lines with the same identifiers
       * (i.e. datetime, source, level), this function returns a list of lines suitable for the
       * caller to pass to the logger one after the other.
       */
      virtual vector<string> summary() const;
      /**
       * Prints this Status object to the logger output stream.
       *
       * @see [easylogging++ documentation](https://github.com/easylogging/easyloggingpp#logging-your-own-class)
       */
      virtual void log(el::base::type::ostream_t& os) const;
      //! @}
  };
}  // ::pfasst


#include "pfasst/controller/status_impl.hpp"

#endif  // _PFASST__CONTROLLER__STATUS_HPP_
