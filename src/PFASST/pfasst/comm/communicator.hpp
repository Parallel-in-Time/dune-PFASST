#ifndef _PFASST__COMM__COMMUNICATOR_HPP_
#define _PFASST__COMM__COMMUNICATOR_HPP_

#include <memory>

#include "pfasst/controller/status.hpp"


namespace pfasst
{
  /**
   * @defgroup Communicators Communicators
   *   Communicators are required for the parallel execution of _PFASST_ and provide ways of passing
   *   data between multiple processors.
   * @ingroup Assistance
   */
  namespace comm
  {
    /**
     * @ingroup Communicators
     */
    class Communicator
      : public std::enable_shared_from_this<Communicator>
    {
      public:
        Communicator() = default;
        Communicator(const Communicator& other) = default;
        Communicator(Communicator&& other) = default;
        virtual ~Communicator() = default;
        Communicator& operator=(const Communicator& other) = default;
        Communicator& operator=(Communicator&& other) = default;

        virtual size_t get_size() const;
        virtual size_t get_rank() const;
        virtual size_t get_root() const;

        virtual bool is_first() const;
        virtual bool is_last() const;

        virtual void cleanup(const bool discard = false);
        virtual void abort(const int& err_code);

        virtual bool probe(const int src_rank, const int tag);

        virtual void send(const double* const data, const int count, const int dest_rank, const int tag);
        virtual void send_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag);

        virtual void isend(const double* const data, const int count, const int dest_rank, const int tag);
        virtual void isend_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag);

        virtual void recv(double* data, const int count, const int src_rank, const int tag);
        virtual void recv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag);

        virtual void irecv(double* data, const int count, const int src_rank, const int tag);
        virtual void irecv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag);

        virtual void bcast(double* data, const int count, const int root_rank);
    };
  }  // ::pfasst::comm
}  // ::pfasst

#include "pfasst/comm/communicator_impl.hpp"

#endif  // _PFASST__COMM__COMMUNICATOR_HPP_
