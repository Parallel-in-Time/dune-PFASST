#include "pfasst/encap/eigen3_vector.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <type_traits>
using std::shared_ptr;

#include "pfasst/logging.hpp"
#include "pfasst/util.hpp"


namespace pfasst
{
  namespace encap
  {
    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(const size_t size)
      : _data(size)
    {
      this->zero();
    }

    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(const typename EncapsulationTrait::data_t& data)
      : Encapsulation<EncapsulationTrait>(data.rows())
    {
      this->data() = data;
    }

    template<class EncapsulationTrait>
    Encapsulation<EncapsulationTrait>&
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::operator=(const typename EncapsulationTrait::data_t& data)
    {
      this->data() = data;
      return *this;
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::data_t&
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::data()
    {
      return this->_data;
    }

    template<class EncapsulationTrait>
    const typename EncapsulationTrait::data_t&
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_data() const
    {
      return this->_data;
    }

    template<class EncapsulationTrait>
    size_t
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_total_num_dofs() const
    {
      assert(this->get_data().rows() == this->get_data().cols());
      return this->get_data().rows() * this->get_data().cols();
    }

    template<class EncapsulationTrait>
    std::array<int, EncapsulationTrait::DIM>
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_dimwise_num_dofs() const
    {
      std::array<int, EncapsulationTrait::DIM> dimwise_ndofs;
      switch (EncapsulationTrait::DIM) {
        case 1:
          dimwise_ndofs.fill((int)this->get_total_num_dofs());
          break;
        case 2:
          assert(this->get_data().cols() == this->get_data().rows());
          dimwise_ndofs.fill((int)this->get_data().cols());
          break;
        default:
          ML_CLOG(FATAL, "ENCAP", "unsupported spatial dimension: " << EncapsulationTrait::DIM);
          throw std::runtime_error("unsupported spatial dimension");
      }

      return dimwise_ndofs;
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::zero()
    {
      this->data().setZero();
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::scaled_add(const typename EncapsulationTrait::time_t& a,
                                   const shared_ptr<Encapsulation<EncapsulationTrait>> y)
    {
      this->data() += a * y->get_data();
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::spatial_t
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::norm0() const
    {
      return this->get_data().template lpNorm<Eigen::Infinity>();
    }

    template<class EncapsulationTrait>
    template<class CommT>
    bool
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::probe(shared_ptr<CommT> comm, const int src_rank, const int tag)
    {
      return comm->probe(src_rank, tag);
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::send(shared_ptr<CommT> comm, const int dest_rank,
                              const int tag, const bool blocking)
    {
      ML_CVLOG(2, "ENCAP", "sending data: " << this->get_data());
      if (blocking) {
        comm->send(this->get_data().data(), this->get_data().cols() * this->get_data().rows(), dest_rank, tag);
      } else {
        comm->isend(this->get_data().data(), this->get_data().cols() * this->get_data().rows(), dest_rank, tag);
      }
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::recv(shared_ptr<CommT> comm, const int src_rank,
                              const int tag, const bool blocking)
    {
      if (blocking) {
        comm->recv(this->data().data(), this->get_data().cols() * this->get_data().rows(), src_rank, tag);
      } else {
        comm->irecv(this->data().data(), this->get_data().cols() * this->get_data().rows(), src_rank, tag);
      }
      ML_CVLOG(2, "ENCAP", "received data: " << this->get_data());
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::bcast(shared_ptr<CommT> comm, const int root_rank)
    {
      comm->bcast(this->data().data(), this->get_data().cols() * this->get_data().rows(), root_rank);
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait, 
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::log(el::base::type::ostream_t& os) const
    {
      os << "EigenVector(" << this->get_data() << ")";
    }


    template<class EncapsulationTrait>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::EncapsulationFactory(const size_t size)
      : _size(size)
    {}

    template<class EncapsulationTrait>
    shared_ptr<Encapsulation<EncapsulationTrait>>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::create() const
    {
      return std::make_shared<Encapsulation<EncapsulationTrait>>(this->size());
    }

    template<class EncapsulationTrait>
    void
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::set_size(const size_t& size)
    {
      this->_size = size;
    }

    template<class EncapsulationTrait>
    size_t
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<eigen3_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::size() const
    {
      return this->_size;
    }
  }  // ::pfasst::encap
}  // ::pfasst
