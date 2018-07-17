#include "pfasst/encap/encapsulation.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace encap
  {
    template<class EncapsulationTrait>
    shared_ptr<Encapsulation<EncapsulationTrait>>
    axpy(const typename EncapsulationTrait::time_t& a,
         const shared_ptr<Encapsulation<EncapsulationTrait>> x,
         const shared_ptr<Encapsulation<EncapsulationTrait>> y)
    {
      shared_ptr<Encapsulation<EncapsulationTrait>> result = \
        std::make_shared<Encapsulation<EncapsulationTrait>>(*y);
      result->scaled_add(a, x);
      return result;
    }

    template<class EncapsulationTrait>
    void
    mat_apply(vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x,
              const typename EncapsulationTrait::time_t& a,
              const Matrix<typename EncapsulationTrait::time_t>& mat,
              const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& y,
              const bool zero_vec_x)
    {
       int rank, num_pro;
       MPI_Comm_rank(MPI_COMM_WORLD, &rank );
      CLOG_IF(x.size() != (size_t)mat.rows(), WARNING, "ENCAP")
        << "size of result vector (" << x.size()
        << ") does not match result of matrix-vector multiplication (" << mat.rows() << ")";

      if (zero_vec_x) {
        std::for_each(x.begin(), x.end(),
                 [](shared_ptr<Encapsulation<EncapsulationTrait>> xi) {
                   xi->zero();
                });
      }

      const size_t cols = mat.cols();
      const size_t rows = mat.rows();

      for (size_t n = 0; n < rows; ++n) {
        for (size_t m = 0; m < cols; ++m) {
         
          x[n]->scaled_add(a * mat(n, m), y[m]);

        }
      }
      //MPI_Barrier(MPI_COMM_WORLD);//std::exit(0);
    }

    template<class EncapsulationTrait>
    vector<shared_ptr<Encapsulation<EncapsulationTrait>>>
    mat_mul_vec(const typename EncapsulationTrait::time_t& a,
                const Matrix<typename EncapsulationTrait::time_t>& mat,
                const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x)
    {
      assert((size_t)mat.cols() == x.size());
      const size_t rows = mat.rows();

      // initialize result vector of encaps
      vector<shared_ptr<Encapsulation<EncapsulationTrait>>> result(rows);
      for(auto& ri : result) {
        ri = std::make_shared<Encapsulation<EncapsulationTrait>>();
        ri->data() = x[0]->get_data();
        ri->zero();
      }

      mat_apply(result, a, mat, x, false);
      return result;
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::spatial_t
    norm0(const shared_ptr<Encapsulation<EncapsulationTrait>> x)
    {
      return x->norm0();
    }
  }  // ::pfasst::encap
}  // ::pfasst
