#include "pfasst/controller/status.hpp"

#include <cstddef>  // offsetof
#include <memory>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>
using std::shared_ptr;
using std::string;
using std::vector;

#include <pfasst/logging.hpp>


namespace pfasst
{
//#ifdef WITH_MPI
  template<typename precision>
  void
  Status<precision>::create_mpi_datatype()
  {
    const int COUNT = 11;
    int blocks[COUNT] = {
      sizeof(PrimaryState),   // primary_state
      sizeof(SecondaryState), // secondary_state
      sizeof(size_t),         // step
      sizeof(size_t),         // num_steps
      sizeof(size_t),         // iteration
      sizeof(size_t),         // max_iterations
      sizeof(precision),      // time
      sizeof(precision),      // dt
      sizeof(precision),      // t_end
      sizeof(precision),      // abs_res_norm
      sizeof(precision)       // rel_res_norm
    };
    MPI_Aint displ[COUNT] = {
      offsetof(StatusDetail<precision>, primary_state),
      offsetof(StatusDetail<precision>, secondary_state),
      offsetof(StatusDetail<precision>, step),
      offsetof(StatusDetail<precision>, num_steps),
      offsetof(StatusDetail<precision>, iteration),
      offsetof(StatusDetail<precision>, max_iterations),
      offsetof(StatusDetail<precision>, time),
      offsetof(StatusDetail<precision>, dt),
      offsetof(StatusDetail<precision>, t_end),
      offsetof(StatusDetail<precision>, abs_res_norm),
      offsetof(StatusDetail<precision>, rel_res_norm)
    };
    MPI_Datatype types[COUNT] = {
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE,
      MPI_BYTE
    };

    ML_CVLOG(1, "DEFAULT", "creating MPI Data Type for Status");
    int err = MPI_Type_create_struct(COUNT, blocks, displ, types, &(status_data_type));
    assert(err == MPI_SUCCESS);
    err = MPI_Type_commit(&(status_data_type));
    assert(err == MPI_SUCCESS);
  }

  template<typename precision>
  void
  Status<precision>::free_mpi_datatype()
  {
    int err = MPI_Type_free(&(status_data_type));
    assert(err == MPI_SUCCESS);
  }
//#endif

  /**
   * @details It actually replaces the internal data storage by a new object throwing away current
   *   data.
   */
  template<typename precision>
  void
  Status<precision>::clear()
  {
    this->_detail = StatusDetail<precision_t>();
  }

  /**
   * @details The following values are reset to their initial/default values:
   *   * primary state
   *   * iteration
   *   * abs_res_norm
   *   * rel_res_norm
   */
  template<typename precision>
  void
  Status<precision>::reset()
  {
    this->set_primary_state(PrimaryState::UNKNOWN_PRIMARY);
    this->iteration() = 0;
    this->abs_res_norm() = 0.0;
    this->rel_res_norm() = 0.0;
  }

  template<typename precision>
  size_t
  Status<precision>::get_step() const
  {
    return this->_detail.step;
  }

  template<typename precision>
  size_t&
  Status<precision>::step()
  {
    return this->_detail.step;
  }

  template<typename precision>
  size_t
  Status<precision>::get_num_steps() const
  {
    return this->_detail.num_steps;
  }

  template<typename precision>
  size_t&
  Status<precision>::num_steps()
  {
    return this->_detail.num_steps;
  }

  template<typename precision>
  size_t
  Status<precision>::get_iteration() const
  {
    return this->_detail.iteration;
  }

  template<typename precision>
  size_t& Status<precision>::iteration()
  {
    return this->_detail.iteration;
  }

  template<typename precision>
  size_t
  Status<precision>::get_max_iterations() const
  {
    return this->_detail.max_iterations;
  }

  template<typename precision>
  size_t&
  Status<precision>::max_iterations()
  {
    return this->_detail.max_iterations;
  }

  template<typename precision>
  precision
  Status<precision>::get_time() const
  {
    return this->_detail.time;
  }

  template<typename precision>
  precision&
  Status<precision>::time()
  {
    return this->_detail.time;
  }

  template<typename precision>
  precision
  Status<precision>::get_dt() const
  {
    return this->_detail.dt;
  }

  template<typename precision>
  precision&
  Status<precision>::dt()
  {
    return this->_detail.dt;
  }

  template<typename precision>
  precision
  Status<precision>::get_t_end() const
  {
    return this->_detail.t_end;
  }

  template<typename precision>
  precision&
  Status<precision>::t_end()
  {
    return this->_detail.t_end;
  }

  template<typename precision>
  PrimaryState
  Status<precision>::get_primary_state() const
  {
    return this->_detail.primary_state;
  }

  template<typename precision>
  void
  Status<precision>::set_primary_state(const PrimaryState& state)
  {
    this->_detail.primary_state = state;
    this->_detail.secondary_state = SecondaryState::UNKNOWN_SECONDARY;
  }

  template<typename precision>
  SecondaryState
  Status<precision>::get_secondary_state() const
  {
    return this->_detail.secondary_state;
  }

  /**
   * @throws std::runtime_error if @p state can not be combined with the currently stored
   *                            PrimaryState.
   */
  template<typename precision>
  void
  Status<precision>::set_secondary_state(const SecondaryState& state)
  {
    if (validate_state_combination(this->get_primary_state(), state)) {
      this->_detail.secondary_state = state;
    } else {
      ML_CLOG(FATAL, "DEFAULT", "Invalid combination of primary ("
                                << (+this->get_primary_state())._to_string()
                                << ") and secondary state ("
                                << (+state)._to_string() << ")");
      throw std::runtime_error("invalid combination of primary and secondary state");
    }
  }

  template<typename precision>
  precision
  Status<precision>::get_abs_res_norm() const
  {
    return this->_detail.abs_res_norm;
  }

  template<typename precision>
  precision&
  Status<precision>::abs_res_norm()
  {
    return this->_detail.abs_res_norm;
  }

  template<typename precision>
  precision
  Status<precision>::get_rel_res_norm() const
  {
    return this->_detail.rel_res_norm;
  }

  template<typename precision>
  precision&
  Status<precision>::rel_res_norm()
  {
    return this->_detail.rel_res_norm;
  }

  template<typename precision>
  template<typename CommT>
  bool
  Status<precision>::probe(shared_ptr<CommT> comm, const int src_rank, const int tag)
  {
    return comm->probe(src_rank, tag);
  }

  template<typename precision>
  template<typename CommT>
  void
  Status<precision>::send(shared_ptr<CommT> comm, const int dest_rank, const int tag,
                          const bool blocking)
  {
    if (blocking) {
      comm->send_status(&(this->_detail), 1, dest_rank, tag);
    } else {
      comm->isend_status(&(this->_detail), 1, dest_rank, tag);
    }
  }

  template<typename precision>
  template<typename CommT>
  void
  Status<precision>::recv(shared_ptr<CommT> comm, const int src_rank, const int tag,
                          const bool blocking)
  {
    if (blocking) {
      comm->recv_status(&(this->_detail), 1, src_rank, tag);
    } else {
      comm->irecv_status(&(this->_detail), 1, src_rank, tag);
    }
  }

  template<typename precision>
  vector<string>
  Status<precision>::summary() const
  {
    vector<string> out;
    out.push_back("Number Iterations: " + std::to_string(this->get_iteration()));
    {
      std::stringstream os;
      //os << "Absolute Residual: " << LOG_FLOAT << this->get_abs_res_norm();
      os << "Absolute Residual: " << this->get_abs_res_norm();
      out.push_back(os.str());
    }
    {
      std::stringstream os;
      //os << "Relative Residual: " << LOG_FLOAT << this->get_rel_res_norm();
      os << "Absolute Residual: " << this->get_abs_res_norm();
      out.push_back(os.str());
    }
    return out;
  }

  template<typename precision>
  void
  Status<precision>::log(el::base::type::ostream_t& os) const
  {
    os << "Status(" << this->_detail << ")";
  }
}  // ::pfasst
