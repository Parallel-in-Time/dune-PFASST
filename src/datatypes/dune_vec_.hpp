#ifndef _PFASST_DUNE_VECTOR_HPP_
#define _PFASST_DUNE_VECTOR_HPP_

#include <memory>
//#include <vector>
using namespace std;

#ifdef WITH_MPI
#include "pfasst/mpi_communicator.hpp"
using namespace pfasst::mpi;
#endif

#include "pfasst/encap/encapsulation.hpp"


#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

using Dune::BlockVector;
using Dune::FieldVector;

namespace pfasst
{
  namespace encap
  {
    /**
     * @tparam scalar
     *     precision and numerical type of the data values
     * @tparam time
     *     precision of the time points; defaults to pfasst::time_precision
     */
    template<typename scalar, typename time = time_precision>
    class Dune_VectorEncapsulation
      : public BlockVector<FieldVector<scalar,1>>, //public BlockVector<FieldVector<SpatialPrecision,1>>>
      //: public vector<scalar>,
        public Encapsulation<time>
    {
      public:

        //! @{
        Dune_VectorEncapsulation(const size_t size);

        /**
         * Copy constuctor.
         *
         * @note delegated to sdt::vector<scalar>
         */
        Dune_VectorEncapsulation(const Dune_VectorEncapsulation<scalar, time>& other);

        /**
         * @throws std::bad_cast
         *     if `other` can not be transformed into pfasst::encap::Dune_VectorEncapsulation via
         *     `dynamic_cast`
         */
        Dune_VectorEncapsulation(const Encapsulation<time>& other);

        /**
         * Move constructor.
         *
         * @note delegated to std::vector<scalar>
         */
        Dune_VectorEncapsulation(Dune_VectorEncapsulation<scalar, time>&& other);

        /**
         * @throws std::bad_cast
         *     if `other` can not be transformed into pfasst::encap::Dune_VectorEncapsulation via
         *     `dynamic_cast`
         */
        Dune_VectorEncapsulation(Encapsulation<time>&& other);

        virtual ~Dune_VectorEncapsulation();
        //! @}

        //! @{
        virtual void zero() override;
        virtual void copy(shared_ptr<const Encapsulation<time>> x) override;
        virtual void copy(shared_ptr<const Dune_VectorEncapsulation<scalar, time>> x);
        //! @}

        //! @{
        virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override;
        virtual void saxpy(time a, shared_ptr<const Dune_VectorEncapsulation<scalar, time>> x);

        /**
         * @note In case any of the elements of `dst` or `src` can not be transformed via
         *     `dynamic_cast` into pfasst::encap::Dune_VectorEncapsulation std::abort is called.
         */
        virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst,
                               time a, Matrix<time> mat,
                               vector<shared_ptr<Encapsulation<time>>> src,
                               bool zero = true) override;
        virtual void mat_apply(vector<shared_ptr<Dune_VectorEncapsulation<scalar, time>>> dst,
                               time a, Matrix<time> mat,
                               vector<shared_ptr<Dune_VectorEncapsulation<scalar, time>>> src,
                               bool zero = true);

        /**
         * Maximum norm of contained elements.
         *
         * This uses std::max with custom comparison function.
         */
        virtual time norm0() const override;
        //! @}

#ifdef WITH_MPI
        //! @{
        MPI_Request recv_request = MPI_REQUEST_NULL;
        MPI_Request send_request = MPI_REQUEST_NULL;
        //! @}

        //! @{
        inline MPICommunicator& as_mpi(ICommunicator* comm)
        {
          auto mpi = dynamic_cast<MPICommunicator*>(comm);
          assert(mpi);
          return *mpi;
        }
        //! @}

        //! @{
        virtual void post(ICommunicator* comm, int tag) override;
        virtual void recv(ICommunicator* comm, int tag, bool blocking) override;
        virtual void send(ICommunicator* comm, int tag, bool blocking) override;
        virtual void broadcast(ICommunicator* comm) override;
        //! @}
#endif

    };

    /**
     * @tparam scalar
     *     precision and numerical type of the data values
     * @tparam time
     *     precision of the time points; defaults to pfasst::time_precision
     */
    template<typename scalar, typename time = time_precision>
    class Dune_VectorFactory
      : public EncapFactory<time>
    {
      protected:
        size_t size;

      public:
        Dune_VectorFactory(const size_t size);
        virtual shared_ptr<Encapsulation<time>> create(const EncapType) override;
        size_t dofs() const;
    };

    template<typename scalar, typename time = time_precision>
    Dune_VectorEncapsulation<scalar,time>& as_vector(shared_ptr<Encapsulation<time>> x);

    template<typename scalar, typename time = time_precision>
    const Dune_VectorEncapsulation<scalar,time>& as_vector(shared_ptr<const Encapsulation<time>> x);
  }  // ::pfasst::encap
}  // ::pfasst

#include "dune_vec_impl.hpp"

#endif