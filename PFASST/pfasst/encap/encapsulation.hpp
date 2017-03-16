#ifndef _PFASST__ENCAP__INTERFACE_HPP_
#define _PFASST__ENCAP__INTERFACE_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/encap/traits.hpp"


namespace pfasst
{
  /**
   * @defgroup Encapsulation Encapsulation
   *   Encapsulations provide a uniform wrapper around user data for the Controllers and Sweepers to
   *   work on.
   * @ingroup Assistance
   */

  /**
   * Dealing with user data.
   *
   * @ingroup Encapsulation
   */
  namespace encap
  {
    template<
      class EncapsulationTrait,
      typename Enabled = void
    >
    class Encapsulation;

    template<
      class EncapsulationTrait,
      typename Enabled = void
    >
    class EncapsulationFactory;


    /**
     * Computes \\( a * x + y \\) .
     *
     * @tparam    EncapsulationTrait type traits for Encapsulations @p x and @p y
     * @param[in] a                  scalar factor to scale @p x with
     * @param[in] x                  first Encapsulation to be scaled by @p a
     * @param[in] y                  second Encapsulation to be added to \\( a*x \\)
     * @returns result of \\( a * x + y \\)
     *
     * @see pfasst::encap_traits
     *  for instances of @p EncapsulationTrait
     *
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait
    >
    shared_ptr<Encapsulation<EncapsulationTrait>>
    axpy(const typename EncapsulationTrait::time_type& a,
         const shared_ptr<Encapsulation<EncapsulationTrait>> x,
         const shared_ptr<Encapsulation<EncapsulationTrait>> y);

    /**
     * Computes scaled matrix-vector product \\( x += aMy \\) .
     *
     * The matrix-vector product of @p mat and @p y is scaled by @p a and added onto @p x.
     *
     * @tparam        EncapsulationTrait type traits for Encapsulations @p x and @p y
     * @param[in,out] x                  target Encapsulation
     * @param[in]     a                  scalar factor to scale the matrix-vector product of @p mat and @p y with
     * @param[in]     matrix             matrix
     * @param[in]     y                  second Encapsulation used as vector in \\( aMy \\)
     * @param[in]     zero_vec_x         if `true` @p x is zeroed before adding the result of the scaled matrix-vector
     *                                   product
     *
     * @see pfasst::encap_traits
     *  for instances of @p EncapsulationTrait
     * @see pfasst::encap::mat_mul_vec()
     *  for a version not modifying any input data
     *
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait
    >
    void
    mat_apply(vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x,
              const typename EncapsulationTrait::time_type& a,
              const Matrix<typename EncapsulationTrait::time_type>& matrix,
              const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& y,
              const bool zero_vec_x = true);

    /**
     * Computes scaled matrix-vector product \\( aMx \\).
     *
     * @tparam    EncapsulationTrait type traits for Encapsulations @p x
     * @param[in] x                  data used as vector
     * @param[in] a                  scalar factor to scale the matrix-vector product of @p mat and @p x with
     * @param[in] matrix             matrix
     * @returns result of \\( aMx \\)
     *
     * @see pfasst::encap_traits
     *  for instances of @p EncapsulationTrait
     * @see pfasst::encap::mat_apply()
     *  for an in-place version
     *
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait
    >
    vector<shared_ptr<Encapsulation<EncapsulationTrait>>>
    mat_mul_vec(const typename EncapsulationTrait::time_type& a,
                const Matrix<typename EncapsulationTrait::time_type>& matrix,
                const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x);

    /**
     * Computes maximums norm of @p x.
     *
     * Computes the maximums or infinity norm of @p x, e.g. \\( \\|x\\|_{\\inf} = max(|x_0|, \\dots, |x_n|) \\).
     *
     * @tparam    EncapsulationTrait type traits for Encapsulations @p x
     * @param[in] x                  data vector to compute infinity/maximums norm of
     * @returns maximums norm
     *
     * @see Encapsulation::norm0()
     *
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait
    >
    typename EncapsulationTrait::spatial_type
    norm0(const shared_ptr<Encapsulation<EncapsulationTrait>> x);


    /**
     * Encapsulations are the way _PFASST_ can handle arbitrary user data.
     *
     * @tparam EncapsulationTrait type trait describing encapsulated data
     * @tparam Enabled            utility type for template specializations
     *
     * @see pfasst::encap::encap_traits
     *  for instances of @p EncapsulationTrait
     *
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait,
      typename Enabled
    >
    class Encapsulation
      :   public std::enable_shared_from_this<Encapsulation<EncapsulationTrait, Enabled>>
        , el::Loggable
    {
      public:
        using traits = EncapsulationTrait;
        using factory_t = EncapsulationFactory<traits>;

      static_assert(std::is_arithmetic<typename traits::time_t>::value,
                    "time precision must be an arithmetic type");
      static_assert(std::is_arithmetic<typename traits::spatial_t>::value,
                    "spatial precision must be an arithmetic type");
      static_assert(std::is_constructible<typename traits::data_t>::value,
                    "Data Type must be constructible");
      static_assert(std::is_default_constructible<typename traits::data_t>::value,
                    "Data Type must be default constructible");
      static_assert(std::is_destructible<typename traits::data_t>::value,
                    "Data Type must be destructible");
      static_assert(std::is_assignable<typename traits::data_t, typename traits::data_t>::value,
                    "Data Type must be assignable");

      //! @cond static_warnings
      STATIC_WARNING(std::is_move_constructible<typename traits::data_t>::value,
                     "Data Type should be move constructible");
      STATIC_WARNING(std::is_copy_constructible<typename traits::data_t>::value,
                     "Data Type should be copy constructible");
      STATIC_WARNING(std::is_move_assignable<typename traits::data_t>::value,
                     "Data Type should be move assignable");
      STATIC_WARNING(std::is_copy_assignable<typename traits::data_t>::value,
                     "Data Type should be copy assignable");
      //! @endcond

      protected:
        //! @{
        //! actual storage of encapsulated data
        typename traits::data_t _data;
        //! @}

      public:
        //! @{
        Encapsulation() = default;
        /**
         * Encapsulating existing data by copying.
         *
         * A copy of @p data will be encapsulated in the newly constructed Encapsulation.
         *
         * @param[in] data
         */
        Encapsulation(const typename EncapsulationTrait::data_t& data);
        Encapsulation(const Encapsulation<EncapsulationTrait, Enabled>& other)= default;
        Encapsulation(Encapsulation<EncapsulationTrait, Enabled>&& other) = default;
        virtual ~Encapsulation() = default;
        Encapsulation<EncapsulationTrait, Enabled>& operator=(const typename EncapsulationTrait::data_t& data);
        Encapsulation<EncapsulationTrait, Enabled>& operator=(const Encapsulation<EncapsulationTrait, Enabled>& other)= default;
        Encapsulation<EncapsulationTrait, Enabled>& operator=(Encapsulation<EncapsulationTrait, Enabled>&& other)= default;
        //! @}

        //! @name Accessor
        //! @{
        /**
         * Accessor for encapsulated data for modification.
         *
         * @returns encapsulated data for modification, e.g. as an lvalue
         */
        virtual       typename EncapsulationTrait::data_t& data();
        //! Read-only version of `data()`.
        virtual const typename EncapsulationTrait::data_t& get_data() const;

        /**
         * Total number of spatial degrees of freedom.
         */
        virtual size_t get_total_num_dofs() const;
        /**
         * Total number of spatial degrees of freedom for each dimension.
         *
         * @returns array of DOFs per dimension; length of array matches `traits::DIM`.
         */
        virtual std::array<int, EncapsulationTrait::DIM> get_dimwise_num_dofs() const;

        /**
         * Computes maximums norm of encapsulated data.
         *
         * @returns maximums norm of underlying data
         *
         * @note The implementation is strongly dependent on the encapsulated data type.
         */
        virtual typename EncapsulationTrait::spatial_t norm0() const;

        /**
         * Streams string representation of Encapsulation.
         *
         * @param[in,out] os stream used for output
         *
         * @see [Documentation of easylogging++](https://github.com/easylogging/easyloggingpp#logging-your-own-class)
         *  for details on where this comes from
         */
        virtual void log(el::base::type::ostream_t& os) const override;
        //! @}

        //! @name Modification
        //! @{
        /**
         * Zeros the underlying data without changing its reserved and occupied space.
         */
        virtual void zero();

        /**
         * Adds @p y scaled by @p a onto this Encapsulation.
         *
         * @param[in] a scaling factor
         * @param[in] y other data to be scaled and added onto this one
         */
        virtual void scaled_add(const typename EncapsulationTrait::time_t& a,
                                const shared_ptr<Encapsulation<EncapsulationTrait>> y);
        //! @}

        //! @name Communication
        //! @{
        /**
         * probing for incomming encapsulated data over communicator.
         *
         * @param[in] comm     Communicator used for sending
         * @param[in] src_rank source processor of the data
         * @param[in] tag      accociation of the data
         * @returns `true` if data is available for receiving, `false` otherwise
         */
        template<class CommT>
        bool probe(shared_ptr<CommT> comm, const int src_rank, const int tag);

        /**
         * Sending encapsulated data over communicator.
         *
         * @param[in] comm      Communicator used for sending
         * @param[in] dest_rank target processor of the data
         * @param[in] tag       accociation of the data
         * @param[in] blocking  `true` for blocking sending, `false` for non-blocking communication
         */
        template<class CommT>
        void send(shared_ptr<CommT> comm, const int dest_rank, const int tag, const bool blocking);

        /**
         * Receiving encapsulated data over communicator.
         *
         * @param[in] comm      Communicator used for sending
         * @param[in] src_rank  source processor of the data
         * @param[in] tag       accociation of the data
         * @param[in] blocking  `true` for blocking receiving, `false` for non-blocking communication
         */
        template<class CommT>
        void recv(shared_ptr<CommT> comm, const int src_rank, const int tag, const bool blocking);

        /**
         * Sending encapsulated data over communicator.
         *
         * @param[in] comm      Communicator used for broadcasting
         * @param[in] root_rank source processor of the data
         */
        template<class CommT>
        void bcast(shared_ptr<CommT> comm, const int root_rank);
        //! @}
    };


    /**
     * Utility to create specific Encapsulations with predefined defaults.
     *
     * @tparam EncapsulationTrait type trait describing encapsulated data
     * @tparam Enabled            utility type for template specializations
     *
     * @note Specializations for certain encapsulated data types may define additional member functions for setting
     *  and accessing default values for instantiated Encapsulations.
     *
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait,
      class Enabled
    >
    class EncapsulationFactory
    {
      public:
        using encap_t = Encapsulation<EncapsulationTrait>;

        //! @{
        EncapsulationFactory();
        EncapsulationFactory(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory(EncapsulationFactory<EncapsulationTrait>&& other);
        virtual ~EncapsulationFactory() = default;
        EncapsulationFactory<EncapsulationTrait>& operator=(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory<EncapsulationTrait>& operator=(EncapsulationFactory<EncapsulationTrait>&& other);
        //! @}

        //! @{
        /**
         * Instantiating a new Encapsulation.
         *
         * @returns a newly instantiated Encapsulation with implementation dependent defaults.
         */
        shared_ptr<Encapsulation<EncapsulationTrait>> create() const;
        //! @}
    };
  }  // ::encap
}  // ::pfasst


#include "pfasst/encap/encapsulation_impl.hpp"

#endif  // _PFASST__ENCAP__INTERFACE_HPP_
