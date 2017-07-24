//#include "sp.hpp"

#include <cassert>
#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"



namespace pfasst
{
  namespace contrib
  {
  
    
 

    
    
    
    template<class TransferTraits>
    void
    SpectralTransfer<TransferTraits>::create(std::shared_ptr<std::vector<MatrixType*>> transfer)
    {
	
	  std::shared_ptr<std::vector<MatrixType*>> vecvec(transfer);
          std::cout << "tranfer create " << std::endl;
	  set_matrix(*vecvec->at(0), *vecvec->at(0));
	    
    }    
    
    
    
    

    template<class TransferTraits>
    void
    SpectralTransfer<
            TransferTraits>::set_matrix(Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate, Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict)
    {
        
        std::cout << "set matrix " << std::endl;
	    interpolate_matrix = interpolate;
        
        
        /*for (int i=0; i< interpolate_matrix.N(); i++){
	      for (int j=0; j< interpolate_matrix.M(); j++){
            if(interpolate_matrix.exists(i,j)){	
                std::cout<< interpolate_matrix[i][j] << " " ;
            }else{
                std::cout<< 0 << " " ;
            }

	      }
	      std::cout<< ""<< std::endl;
        }*/
    
	    
	    restrict_matrix   = restrict;

	    for (int i=0; i< restrict_matrix.N(); i++){
	      for (int j=0; j< restrict_matrix.M(); j++){
		if(restrict_matrix.exists(i,j)){	
		  if (restrict_matrix[i][j]==0.5 ) restrict_matrix[i][j]=0;
		}

	      }
	    }
	    
	           /* for (int i=0; i< interpolate_matrix.N(); i++){
	      for (int j=0; j< interpolate_matrix.M(); j++){
            if(interpolate_matrix.exists(i,j)){	
                std::cout<< restrict_matrix[i][j] << " " ;
            }else{
                std::cout<< 0 << " " ;
            }

	      }
	      std::cout<< ""<< std::endl;
        }
        std::exit(0);*/
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
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        fine->data() = coarse->get_data();

      } else {

	
	/*std::cout <<  "interpolate grob" <<  std::endl;
        for (int i=0; i< coarse->data().size(); i++){
          std::cout <<  coarse->data()[i] <<  std::endl;
        }
        std::cout <<  "interpolate " <<  std::endl;*/


        interpolate_matrix.mv(coarse->data(), fine->data());
        //Transfer_matrix.mv(coarse->data(), fine->data());
	/*std::cout <<  "interpolate fein" <<  std::endl;
        for (int i=0; i< fine->data().size(); i++){
          std::cout <<  fine->data()[i] <<  std::endl;
        }
        std::cout <<  "interpolate ende" <<  std::endl;*/
        


        
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
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        coarse->data() = fine->get_data();

      } else {


	/*std::cout <<  "restriktion fein" <<  std::endl;
        for (int i=0; i< fine->data().size(); i++){
          std::cout <<  fine->data()[i] <<  std::endl;
        }
        std::cout <<  "restriction " <<  std::endl;*/
	restrict_matrix.mtv(fine->data(), coarse->data());
    //interpolate_matrix.mtv(fine->data(), coarse->data());
        //Transfer_matrix2.mtv(fine->data(), coarse->data());
        //coarse->data() *= 0.5;
	/*std::cout <<  "restriction grob" <<  std::endl;
        for (int i=0; i< coarse->data().size(); i++){
          std::cout <<  coarse->data()[i] <<  std::endl;
        }
        std::cout <<  "restriction ende " <<  std::endl;*/


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
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        coarse->data() = fine->get_data();

      } else {


	/*std::cout <<  "restriktion fein" <<  std::endl;
        for (int i=0; i< fine->data().size(); i++){
          std::cout <<  fine->data()[i] <<  std::endl;
        }
        std::cout <<  "restriction " <<  std::endl;*/
	interpolate_matrix.mtv(fine->data(), coarse->data());
        //Transfer_matrix2.mtv(fine->data(), coarse->data());
        //coarse->data() *= 0.5;
	/*std::cout <<  "restriction grob" <<  std::endl;
        for (int i=0; i< coarse->data().size(); i++){
          std::cout <<  coarse->data()[i] <<  std::endl;
        }
        std::cout <<  "restriction ende " <<  std::endl;*/


      }
    }
    
    
    
    
    
    
  }  // ::pfasst::contrib
}  // ::pfasst



