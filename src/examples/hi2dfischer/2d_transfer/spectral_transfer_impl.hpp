//#include "sp.hpp"

#include <cassert>
#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"

#include <dune/istl/matrixmatrix.hh>

namespace pfasst
{
  namespace contrib
  {
    
    template<class TransferTraits>
    void
    SpectralTransfer<TransferTraits>::create(std::shared_ptr<std::vector<MatrixType*>> transfer)
    {
	
	  //std::shared_ptr<std::vector<MatrixType*>> vecvec(transfer);
	  //set_matrix(*vecvec->at(0), *vecvec->at(0));
	  set_matrix(*transfer->at(0), *transfer->at(0));
	    
    }    
    
    template<class TransferTraits>
    void
    SpectralTransfer<TransferTraits>::create(std::shared_ptr<fe_manager> FinEl)
    {
	
	  std::shared_ptr<std::vector<MatrixType*>> vecvec(FinEl->get_transfer());
          //std::cout << "tranfer create " << std::endl;
          set_matrix(*vecvec->at(0), *vecvec->at(0));

          
	    
    }     
    
    

    template<class TransferTraits>
    void
    SpectralTransfer<
            TransferTraits>::set_matrix(Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate, Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict)
    {
        
	    interpolate_matrix = interpolate;
        
	    restrict_matrix   = restrict;


	    for (int i=0; i< restrict_matrix.N(); i++){
	      for (int j=0; j< restrict_matrix.M(); j++){
		if(restrict_matrix.exists(i,j)){	
		  if (restrict_matrix[i][j]!= 1) restrict_matrix[i][j]=0;//std::cout << restrict_matrix[i][j]<< " ";//
		  //if (restrict_matrix[i][j]==0.5 ) restrict_matrix[i][j]=0;//std::cout << restrict_matrix[i][j]<< " ";//
		}

	      }
	      //std::cout <<  std::endl;
	    }
	  //std::exit(0);
    }


    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                          shared_ptr<typename TransferTraits::fine_encap_t> fine)
    {
      ML_CVLOG(1, "TRANS", "interpolate data");

      const size_t coarse_ndofs = coarse->get_data().size();
      const size_t fine_ndofs = fine->get_data().size();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        std::cout << "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT" << std::endl;
        fine->data() = coarse->get_data();

      } else {

        interpolate_matrix.mv(coarse->data(), fine->data());
        
      }
    }

    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
    {
      ML_CVLOG(1, "TRANS", "restrict data");

      const size_t coarse_ndofs = coarse->get_data().size();
      const size_t fine_ndofs = fine->get_data().size();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        std::cout << "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT" << std::endl;
        coarse->data() = fine->get_data();

      } else {

	restrict_matrix.mtv(fine->data(), coarse->data());
	coarse->data() *= 0.25;
      }
    }
    
    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits>::restrict_u(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
    {
      ML_CVLOG(1, "TRANS", "restrict data");

      const size_t coarse_ndofs = coarse->get_data().size();
      const size_t fine_ndofs = fine->get_data().size();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        std::cout << "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT" << std::endl;
        coarse->data() = fine->get_data();

      } else {

	interpolate_matrix.mtv(fine->data(), coarse->data());


      }
    }

    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits>::restrict_dune_matrix(Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>> f, Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>> c){

	Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>> neu;
	Dune::matMultMat(neu, f, interpolate_matrix); 
	Dune::transposeMatMultMat(c, interpolate_matrix, neu);

   }       
    
  }  // ::pfasst::contrib
}  // ::pfasst



