#ifndef _PFASST__ENCAP__DUNE_VEC_MULTI_HPP_
#define _PFASST__ENCAP__DUNE_VEC_MULTI_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <leathers/push>
#include <leathers/all>
#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>
#include <leathers/pop>
using Dune::BlockVector;
using Dune::FieldVector;
#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/encap/encapsulation.hpp"


//int nr_of_comp=2;

namespace pfasst
{
  namespace encap
  {


      struct dune_encap_tag
              : public encap_data_tag
      {};

      /**
      * Spatialized Type Traits for encapsulation of std::vector.
      *
      * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
      * @tparam SpatialPrecision the spatial data precision
      *
      * @ingroup Traits
      */
      template<
              class TimePrecision,
              class SpatialPrecision,
              size_t Dim,
              size_t Nr_of_comp
      >
      struct dune_vec_encap_traits
              : public encap_traits<TimePrecision, SpatialPrecision, Dim, BlockVector<FieldVector<SpatialPrecision,Nr_of_comp>>>
      {
          using time_t = TimePrecision;
          using spatial_t = SpatialPrecision;
          using data_t = BlockVector<FieldVector<spatial_t,Nr_of_comp>> ; // tu moramo rucno mijenjati 
          using tag_t = dune_encap_tag;
          using dim_t = std::integral_constant<size_t, Dim>;
          static constexpr size_t  DIM = Dim;
	  static constexpr size_t  NR_OF_COMP = Nr_of_comp;
      };




    /**
     * Specialization of Encapsulation for `std::vector`.
     */
    template<
      class EncapsulationTrait
    >
    class Encapsulation<EncapsulationTrait,
                        typename std::enable_if<
                                   std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
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
     */
    template<
      typename time_precision,
      typename spatial_precision,
      size_t Dim,
      size_t Nr_of_comp
    >
    using DuneEncapsulation = Encapsulation<dune_vec_encap_traits<time_precision, spatial_precision, Dim, Nr_of_comp>>;


    template<
      class EncapsulationTrait
    >
    class EncapsulationFactory<EncapsulationTrait,
                               typename std::enable_if<
                                          std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
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
    template<typename T>
    static string join(const BlockVector<T>& vec, const string& sep)
    {
#ifndef PFASST_NO_LOGGING
      std::stringstream os;
      os << "[" << LOG_FLOAT;
      for (size_t i = 0; i < vec.size() - 1; ++i) {
        os << vec[i] << sep;
      }
      os << vec[vec.size()-1];
      os << "]";
      return os.str();
#else
      UNUSED(vec); UNUSED(sep);
    return "";
#endif
    }
}  // ::pfasst


#include "pfasst/encap/dune_vec_multi_impl.hpp"

#endif  // _PFASST__ENCAP__DUNE_VEC_HPP_
