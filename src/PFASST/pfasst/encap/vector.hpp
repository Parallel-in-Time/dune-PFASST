#ifndef _PFASST__ENCAP__VECTOR_HPP_
#define _PFASST__ENCAP__VECTOR_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  namespace encap
  {
    /**
     * Specialization of Encapsulation for `std::vector`.
     *
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait
    >
    class Encapsulation<EncapsulationTrait,
                        typename std::enable_if<
                                   std::is_same<vector_encap_tag, typename EncapsulationTrait::tag_t>::value
                                 >::type>
      :   public std::enable_shared_from_this<Encapsulation<EncapsulationTrait>>
        , public el::Loggable
    {
      public:
        using traits = EncapsulationTrait;
        using factory_t = EncapsulationFactory<traits>;

      protected:
        typename traits::data_t _data;

      public:
        explicit Encapsulation(const size_t size = 0);
        Encapsulation(const typename EncapsulationTrait::data_t& data);
        Encapsulation<EncapsulationTrait>& operator=(const typename EncapsulationTrait::data_t& data);

        virtual       typename EncapsulationTrait::data_t& data();
        virtual const typename EncapsulationTrait::data_t& get_data() const;
        virtual size_t get_total_num_dofs() const;
        // assuming square-shaped space
        virtual std::array<int, EncapsulationTrait::DIM> get_dimwise_num_dofs() const;

        virtual void zero();
        virtual void scaled_add(const typename EncapsulationTrait::time_t& a,
                               const shared_ptr<Encapsulation<EncapsulationTrait>> y);

        virtual typename EncapsulationTrait::spatial_t norm0() const;

        template<class CommT>
        bool probe(shared_ptr<CommT> comm, const int src_rank, const int tag);
        template<class CommT>
        void send(shared_ptr<CommT> comm, const int dest_rank, const int tag, const bool blocking);
        template<class CommT>
        void recv(shared_ptr<CommT> comm, const int src_rank, const int tag, const bool blocking);
        template<class CommT>
        void bcast(shared_ptr<CommT> comm, const int root_rank);

        virtual void log(el::base::type::ostream_t& os) const override;
    };

    /**
     * Shortcut for encapsulation of `std::vector` data types.
     *
     * @ingroup Encapsulation
     */
    template<
      typename time_precision,
      typename spatial_precision,
      size_t Dim
    >
    using VectorEncapsulation = Encapsulation<vector_encap_traits<time_precision, spatial_precision, Dim>>;


    /**
     * @ingroup Encapsulation
     */
    template<
      class EncapsulationTrait
    >
    class EncapsulationFactory<EncapsulationTrait,
                               typename std::enable_if<
                                          std::is_same<vector_encap_tag, typename EncapsulationTrait::tag_t>::value
                                        >::type>
      : public std::enable_shared_from_this<EncapsulationFactory<EncapsulationTrait>>
    {
      protected:
        size_t _size;

      public:
        explicit EncapsulationFactory(const size_t size = 0);
        EncapsulationFactory(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory(EncapsulationFactory<EncapsulationTrait>&& other);
        virtual ~EncapsulationFactory() = default;
        EncapsulationFactory<EncapsulationTrait>& operator=(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory<EncapsulationTrait>& operator=(EncapsulationFactory<EncapsulationTrait>&& other);

        virtual shared_ptr<Encapsulation<EncapsulationTrait>> create() const;

        virtual void set_size(const size_t& size);
        virtual size_t size() const;
    };
  }  // ::pfasst::encap
}  // ::pfasst


#include "pfasst/encap/vector_impl.hpp"

#endif  // _PFASST__ENCAP__VECTOR_HPP_
