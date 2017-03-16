#include "pfasst/comm/mpi_p2p.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

#include "pfasst/logging.hpp"


MAKE_LOGGABLE(MPI_Status, mpi_status, os)
{
  if (   mpi_status.MPI_TAG == MPI_ANY_TAG
         && mpi_status.MPI_SOURCE == MPI_ANY_SOURCE
         && mpi_status.MPI_ERROR == MPI_SUCCESS) {
    os << "MPI_Status(empty)";
  } else {
    char err_str[MPI_MAX_ERROR_STRING];
    int err_len = 0;
    int err = MPI_Error_string(mpi_status.MPI_ERROR, err_str, &err_len);
    pfasst::comm::check_mpi_error(err);
    os << "MPI_Status(source=" << std::to_string(mpi_status.MPI_SOURCE) << ", "
       << "tag=" << std::to_string(mpi_status.MPI_TAG) << ", "
       << "error=" << std::string(err_str, err_len) << ")";
  }
  return os;
}


namespace pfasst
{
  namespace comm
  {
    std::string error_from_code(const int err_code)
    {
      char err_str[MPI_MAX_ERROR_STRING];
      int err_len = 0;
      int err = MPI_Error_string(err_code, err_str, &err_len);
      check_mpi_error(err);
      return std::string(err_str, err_len) + " (code=" + std::to_string(err_code) + ")";
    }


    MPI_Status MPI_Status_factory()
    {
      MPI_Status stat;
      stat.MPI_ERROR = MPI_SUCCESS;
      stat.MPI_SOURCE = MPI_ANY_SOURCE;
      stat.MPI_TAG = MPI_ANY_TAG;
      return stat;
    }

    void check_mpi_error(const int err_code)
    {
      if (err_code != MPI_SUCCESS) {
        std::string err_msg = error_from_code(err_code);
        ML_CLOG(ERROR, "COMM_P2P", "MPI encountered an error: " << err_msg);
        throw std::runtime_error("MPI encountered an error: " + err_msg);
      }
    }


    MpiP2P::MpiP2P(MPI_Comm comm)
      :   _comm(comm)
        , _requests(0)
    {
      log::add_custom_logger("COMM_P2P");

      // get communicator's size and processors rank
      MPI_Comm_size(this->_comm, &(this->_size));
      MPI_Comm_rank(this->_comm, &(this->_rank));

      // get communicator's name (if available)
      int len = -1;
      char buff[MPI_MAX_OBJECT_NAME];
      int err = MPI_Comm_get_name(this->_comm, buff, &len);
      check_mpi_error(err);
      if (len > 0) {
        this->_name = std::string(buff, len);
      }
    }

    MpiP2P::~MpiP2P()
    {
      this->cleanup(true);
    }

    size_t MpiP2P::get_size() const
    {
      assert(this->_size > 0);
      return this->_size;
    }

    size_t MpiP2P::get_rank() const
    {
      assert(this->_rank >= 0);
      return this->_rank;
    }

    std::string MpiP2P::get_name() const
    {
      return this->_name;
    }

    bool MpiP2P::is_first() const
    {
      return (this->get_rank() == this->get_root());
    }

    bool MpiP2P::is_last() const
    {
      return (this->get_rank() == this->get_size() - 1);
    }

    void MpiP2P::cleanup(const bool discard)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "cleaning up " << this->_requests.size() << " dangling request handlers");
      int err = -1;

      for (auto& req : this->_requests) {
        MPI_Status stat = MPI_Status_factory();
        err = MPI_Wait(req.get(), &stat);
        check_mpi_error(err);
        if (*(req.get()) == MPI_REQUEST_NULL) {
          req.reset();
        }
      }

      this->_requests.erase(std::remove(this->_requests.begin(), this->_requests.end(), nullptr),
                            this->_requests.end());

      ML_CLOG(DEBUG, "COMM_P2P", "done");
    }

    void MpiP2P::abort(const int& err_code)
    {
      MPI_Abort(this->_comm, err_code);
    }


    bool MpiP2P::probe(const int src_rank, const int tag)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "probing for incomming message from " << src_rank << " with tag=" << tag);
      MPI_Status stat = MPI_Status_factory();
      int flag = (int)false;
      int err = MPI_Iprobe(src_rank, tag, this->_comm, &flag, &stat);
      check_mpi_error(err);
      ML_CLOG(DEBUG, "COMM_P2P", "probed: " << stat);
      return (bool)flag;
    }


    void MpiP2P::send(const double* const data, const int count, const int dest_rank, const int tag)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "sending " << count << " " << ((count == 1) ? "double" : "doubles")
              << " with tag=" << tag << " to " << dest_rank);

      int err = MPI_Send(mpi_const_cast<void>(data), count, MPI_DOUBLE, dest_rank, tag, this->_comm);
      check_mpi_error(err);
    }

    void MpiP2P::send_status(const StatusDetail<double>* const data, const int count,
                             const int dest_rank, const int tag)
    {
      assert(pfasst::status_data_type != MPI_DATATYPE_NULL);

      ML_CLOG(DEBUG, "COMM_P2P",
              "sending " << count << " " << ((count == 1) ? "Status" : "Stati")
              << " with tag=" << tag << " to " << dest_rank);

      int err = MPI_Send(mpi_const_cast<void>(data), count, status_data_type, dest_rank, tag,
                         this->_comm);
      check_mpi_error(err);
    }


    void MpiP2P::isend(const double* const data, const int count, const int dest_rank, const int tag)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "non-blocking send of " << count << " " << ((count == 1) ? "double" : "doubles")
              << " with tag=" << tag << " to " << dest_rank);

      this->_requests.emplace_back(new MPI_Request(MPI_REQUEST_NULL));

      int err = MPI_Isend(mpi_const_cast<void>(data), count, MPI_DOUBLE, dest_rank, tag,
                          this->_comm, this->_requests.back().get());
      check_mpi_error(err);
    }

    void MpiP2P::isend_status(const StatusDetail<double>* const data, const int count,
                              const int dest_rank, const int tag)
    {
      assert(pfasst::status_data_type != MPI_DATATYPE_NULL);

      ML_CLOG(DEBUG, "COMM_P2P",
              "non-blocking send of " << count << " " << ((count == 1) ? "Status" : "Stati")
              << " with tag=" << tag << " to " << dest_rank);

      this->_requests.emplace_back(new MPI_Request(MPI_REQUEST_NULL));

      int err = MPI_Isend(mpi_const_cast<void>(data), count, status_data_type, dest_rank, tag,
                          this->_comm, this->_requests.back().get());
      check_mpi_error(err);
    }


    void MpiP2P::recv(double* data, const int count, const int dest_rank, const int tag)
    {
      auto stat = MPI_Status_factory();
      ML_CLOG(DEBUG, "COMM_P2P",
              "receiving " << count << " " << ((count == 1) ? "double" : "doubles")
              << " with tag=" << tag << " from " << dest_rank);

      int err = MPI_Recv(data, count, MPI_DOUBLE, dest_rank, tag, this->_comm, &stat);
      check_mpi_error(err);
    }

    void MpiP2P::recv_status(StatusDetail<double>* data, const int count, const int dest_rank, const int tag)
    {
      assert(pfasst::status_data_type != MPI_DATATYPE_NULL);

      auto stat = MPI_Status_factory();
      ML_CLOG(DEBUG, "COMM_P2P",
              "receiving " << count << " " << ((count == 1) ? "Status" : "Stati")
              << " with tag=" << tag << " from " << dest_rank);

      int err = MPI_Recv(data, count, pfasst::status_data_type, dest_rank, tag, this->_comm, &stat);
      check_mpi_error(err);
    }


    void MpiP2P::irecv(double* data, const int count, const int src_rank, const int tag)
    {
//       ML_CLOG(DEBUG, "COMM_P2P",
//               "non-blocking receive of " << count << " " << ((count == 1) ? "double" : "doubles")
//               << " with tag=" << tag << " from " << src_rank);
// 
//       this->_requests.push_back(MPI_REQUEST_NULL);
//       int err = MPI_Irecv(data, count, MPI_DOUBLE, src_rank, tag,
//                           mpi_const_cast<MPI_Comm>(this->_comm), &(this->_requests.back()));
//       check_mpi_error(err);
      throw std::runtime_error("dont irecv");
    }

    void MpiP2P::irecv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag)
    {
//       assert(pfasst::status_data_type != MPI_DATATYPE_NULL);
// 
//       ML_CLOG(DEBUG, "COMM_P2P",
//               "non-blocking receive of " << count << " " << ((count == 1) ? "Status" : "Stati")
//               << " with tag=" << tag << " from " << src_rank);
// 
//       this->_requests.push_back(MPI_REQUEST_NULL);
//       int err = MPI_Irecv(data, count, status_data_type, src_rank, tag,
//                           mpi_const_cast<MPI_Comm>(this->_comm), &(this->_requests.back()));
//       check_mpi_error(err);
      throw std::runtime_error("dont irecv");
    }


    void MpiP2P::bcast(double* data, const int count, const int root_rank)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "broadcasting " << count << " " << ((count == 1) ? "double" : "doubles")
              << " from root " << root_rank);

      int err = MPI_Bcast(data, count, MPI_DOUBLE, root_rank, mpi_const_cast<MPI_Comm>(this->_comm));
      check_mpi_error(err);
    }


  } // ::pfasst::comm
}  // ::pfasst
