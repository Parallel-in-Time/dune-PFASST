/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/spectral_transfer_1d.hpp
 * @since v0.1.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__SPECTRAL_TRANSFER_1D_HPP_
#define _EXAMPLES__ADVEC_DIFF__SPECTRAL_TRANSFER_1D_HPP_

#include <cassert>
#include <cstdlib>
#include <memory>
using namespace std;

#include "../datatypes/dune_vec.hpp"
#include <pfasst/transfer/polynomial.hpp>

#include <dune/fufem/assemblers/basisinterpolationmatrixassembler.hh>

//#include "fe_manager.hpp"
//#include "fft_manager.hpp"
//#include "fftw_workspace_dft1d.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Spectral (FFT) transfer routines.
       *
       * @ingroup AdvectionDiffusion
       */
      template<typename time = pfasst::time_precision>
      class SpectralTransfer
        : public encap::PolyInterpMixin<time>
      {
          using Encapsulation = encap::Encapsulation<double>;

          //FFTManager<FFTWWorkspaceDFT1D<encap::VectorEncapsulation<double>>> _fft;
	  Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate_matrix;
	  Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict_matrix;
	  
        public:
	  
	  void set_matrix(Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate, Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict){
            
	    //std::cout << "**************************************** ***************************************************" << std::endl;    
	    interpolate_matrix = interpolate;
	    
	    restrict_matrix   = restrict;

	    for (int i=0; i< restrict_matrix.N(); i++){
	      for (int j=0; j< restrict_matrix.M(); j++){
		if(restrict_matrix.exists(i,j)){
		  if (restrict_matrix[i][j]==0.5 ) restrict_matrix[i][j]=0;
		}

	      }
	    }
	  }
	  
	  void set_matrix(Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> matrix, int i){
            if(i==0){
	    //std::cout << "**************************************** interpolate ***************************************************" << std::endl;  
	    interpolate_matrix = matrix;
	    }
	    if(i==1){
	    //std::cout << "**************************************** restrict ***************************************************" << std::endl;    
	    restrict_matrix   = matrix;

	    for (int i=0; i< restrict_matrix.N(); i++){
	      for (int j=0; j< restrict_matrix.M(); j++){
		if(restrict_matrix.exists(i,j)){
		  if (restrict_matrix[i][j]==0.5 ) restrict_matrix[i][j]=0;
		}

	      }
	    }
	    }
	  }
	  
	  
	  
	  SpectralTransfer(std::shared_ptr<fe_manager> FinEl, size_t nlevel){
	    
	    //set_matrix(FinEl->get_transfer(0));
	    std::shared_ptr<std::vector<MatrixType*>> vecvec(FinEl->get_transfer());
	    size_t l = FinEl->get_nlevel() - nlevel -1; 
	    //size_t l = nlevel;
	    /*std::cout <<  "transfer konstrukor groesse " << (*vecvec->at(0)).M() <<  std::endl;
	    for (int i=0; i< vecvec->at(0)->N(); i++){
	      for (int j=0; j< (*vecvec->at(0)).M(); j++){
		if(vecvec->at(0)->exists(i,j)){
		  std::cout << ((*vecvec->at(0))[i][j]) << std::endl;
		}
	      }

	    }
	    std::cout <<  "transfer konstrukor" <<  std::endl;*/
	    
	    if(l==0){
	      set_matrix(*vecvec->at(l), *vecvec->at(l));
	    }else if(l==FinEl->get_nlevel() -1){
	      set_matrix(*vecvec->at(l-1), *vecvec->at(l-1));
	    }else{
	      set_matrix(*vecvec->at(l-1), *vecvec->at(l-1)); // l-1, l
	      //erstes l restrict scheint richtig 
	    }  
	    //std::cout <<  "transfer konstrukor" <<  std::endl;
	    
	    /*std::cout <<  "transfer konstrukor" <<  std::endl;
	    Transfer_matrix  = (FinEl->get_transfer(0));
	    Transfer_matrix2 = (FinEl->get_transfer(0));

	    std::cout <<  "transfer konstrukor" <<  std::endl;
	    for (int i=0; i< Transfer_matrix2.N(); i++){
	      for (int j=0; j< Transfer_matrix2.M(); j++){
		if(Transfer_matrix.exists(i,j)){
		  if (Transfer_matrix2[i][j]==0.5 ) Transfer_matrix2[i][j]=0;
		}

	      }
	    }*/
	    
	  }
	  
	  
	  SpectralTransfer(std::shared_ptr<fe_manager> FinEl, size_t nlevel, size_t base_order){
	    typedef Dune::YaspGrid<1> GridType; 
	    typedef DuneFunctionsBasis<Dune::Functions::PQkNodalBasis<GridType::LeafGridView, 1>> B;
	    typedef DuneFunctionsBasis<Dune::Functions::PQkNodalBasis<GridType::LeafGridView, 2>> B2;
	    GridType::LeafGridView gridView = FinEl->get_grid()->leafGridView();
	    B b_coarse(gridView);
	    B2 b_fine(gridView);
	    MatrixType interpol_matrix;
	    assembleBasisInterpolationMatrix(interpol_matrix, b_coarse, b_fine);
	    //std::shared_ptr<std::vector<MatrixType*>> vecvec(FinEl->get_transfer());
	    //size_t l = FinEl->get_nlevel() - nlevel -1; 
	    //set_matrix(*vecvec->at(l-1), *vecvec->at(l-1)); 
	    
	  }
	  
          void interpolate(shared_ptr<Encapsulation> dst,
                           shared_ptr<const Encapsulation> src) override
          {
            auto& fine = encap::as_vector<double, time>(dst);
            auto& coarse = encap::as_vector<double, time>(src);

	   // std::cout << "**************************************** im interpolate ***************************************************" << std::endl;  
	          /*for (int i=0; i< Transfer_matrix.N(); i++){
          for (int j=0; j< Transfer_matrix.M(); j++){
            if(Transfer_matrix.exists(i,j))
            std::cout <<  Transfer_matrix[i][j];
          }
          std::cout <<  "" <<  std::endl;
        }*/
		  
		    /*for (int i=0; i< interpolate_matrix.N(); i++){
          for (int j=0; j< interpolate_matrix.M(); j++){
            if(interpolate_matrix.exists(i,j))
            std::cout <<  interpolate_matrix[i][j];
          }
          std::cout <<  "" <<  std::endl;
        }	  */
		  
		  
	/*std::cout <<  "interpolate grob" <<  std::endl;
        for (int i=0; i< coarse.size(); i++){
          std::cout <<  coarse[i] <<  std::endl;
        }
        std::cout <<  "interpolate ende" <<  std::endl;*/

	//std::cout <<  "interpolate "  << std::endl;	
	//if (interpolate_matrix.M() != coarse.size() || interpolate_matrix.N() != fine.size()) std::cout << "ACHTUNG!!!!" << std::endl;
		//std::cout <<  "interpolate " <<interpolate_matrix.N() << " "<< interpolate_matrix.M()<< " "<< coarse.size() << " " << fine.size() << std::endl;
	    /*for (size_t i = 0; i < fine.N(); i++) {
	      std::cout << "u " << fine[i] << std::endl;
	    }
	    for (int i=0; i< interpolate_matrix.N(); i++){
          for (int j=0; j< interpolate_matrix.M(); j++){
            if(interpolate_matrix.exists(i,j))
            std::cout <<  interpolate_matrix[i][j];
          }
          std::cout <<  "" <<  std::endl;
        }	*/
        interpolate_matrix.mv(coarse, fine);
	//std::cout <<  "interpolate ende" << std::endl;

	
        /*std::cout <<  "interpolate fein" <<  std::endl;
        for (int i=0; i< fine.size(); i++){
          std::cout <<  fine[i] <<  std::endl;
        }
        std::cout <<  "interpolate ende" <<  std::endl;*/
	    

            /*auto* crse_z = this->_fft.get_workspace(crse.size())->forward(crse);
            auto* fine_z = this->_fft.get_workspace(fine.size())->z_ptr();

            for (size_t i = 0; i < fine.size(); i++) {
              fine_z[i] = 0.0;
            }

            double c = 1.0 / crse.size();

            for (size_t i = 0; i < crse.size() / 2; i++) {
              fine_z[i] = c * crse_z[i];
            }

            for (size_t i = 1; i < crse.size() / 2; i++) {
              fine_z[fine.size() - crse.size() / 2 + i] = c * crse_z[crse.size() / 2 + i];
            }

            this->_fft.get_workspace(fine.size())->backward(fine);*/
          }

          void restrict(shared_ptr<Encapsulation> dst,
                        shared_ptr<const Encapsulation> src) override
          {
            auto& fine = encap::as_vector<double, time>(src);
            auto& coarse = encap::as_vector<double, time>(dst);

	    //std::cout << "****************************************im restrict ***************************************************" << std::endl;  
	    /*for (int i=0; i< restrict_matrix.N(); i++){
          for (int j=0; j< restrict_matrix.M(); j++){
            if(restrict_matrix.exists(i,j))
            std::cout <<  restrict_matrix[i][j];
          }
          std::cout <<  "" <<  std::endl;
        }*/
	    
	    
	    /*{std::cout <<  "restrict fein" <<  std::endl;
            for (size_t i = 0; i < fine.size(); i++) {
              std::cout <<  fine[i] <<  std::endl;
            }
            std::cout <<  "restrict ende" <<  std::endl;*/

	    //std::cout <<  "restrict "  << std::endl;
	    //if (restrict_matrix.N() != coarse.size() || restrict_matrix.M() != fine.size()) std::cout << "ACHTUNG!!!!" << std::endl;
	    restrict_matrix.mtv(fine, coarse);
	    //std::cout <<  "restrict " <<restrict_matrix.N() << " "<< restrict_matrix.M()<< " "<< coarse.size() << " " << fine.size() << std::endl;

            //size_t xrat = fine.size() / crse.size();

	    /*std::cout <<  "restrict grob" <<  std::endl;
            for (size_t i = 0; i < coarse.size(); i++) {
              std::cout <<  coarse[i] <<  std::endl;
            }}
            std::cout <<  "restrict ende" <<  std::endl;*/
          }
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#endif // _EXAMPLES__ADVEC_DIFF__SPECTRAL_TRANSFER_1D_HPP_
