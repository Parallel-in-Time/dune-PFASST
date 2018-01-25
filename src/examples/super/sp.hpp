#ifndef _PFASST__TRANSFER__SPECTRAL_TRANSFER_1D_HPP_
#define _PFASST__TRANSFER__SPECTRAL_TRANSFER_1D_HPP_

#include "pfasst/transfer/polynomial.hpp"

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/quadrature.hpp"
#include "pfasst/contrib/fft.hpp"


#include "../../datatypes/pdelab_vec.hpp"

//#include <dune/istl/bcrsmatrix.hh>



const int dim=1;


namespace pfasst
{
  namespace contrib
  {
    /**
     * @ingroup Contributed
     */
    template<
      class TransferTraits
    >
    class SpectralTransfer<TransferTraits> //, typename enable_if<is_same< typename TransferTraits::fine_encap_traits::dim_t, integral_constant<size_t, dim> >::value >::type>
      : public PolynomialTransfer<TransferTraits>
    {
      


      public:
        using traits = TransferTraits;

        /*typedef Dune::YaspGrid<1> Grid;
        typedef Grid::ctype DF;
	typedef Grid::LevelGridView GVl;
	typedef Dune::PDELab::QkLocalFiniteElementMap<GVl,DF,double,1> FEM;
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GVl,FEM,CON,VBE> GFS;*/


	const static int dim=1;
	const int degree =1;
	std::array<int,dim> n;
        Dune::FieldVector<double,dim> h = {1};	      
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        std::shared_ptr<Grid> gridp;
	typedef Grid::LevelGridView GV;
	//typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
	typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF,double, 1> FEM;  
	std::shared_ptr<FEM> fem; 
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	std::shared_ptr<GFS> gfs; 

  	using Zl = Dune::PDELab::Backend::Vector<GFS,double>;
	typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType2;
	typedef Dune::PDELab::NonoverlappingOperator<GFS, MatrixType2, Zl, Zl> NOO;

	std::shared_ptr<NOO> prolong;
	std::shared_ptr<NOO> rest;
 	std::shared_ptr<NOO> inj;



      protected:
        pfasst::contrib::FFT<typename traits::fine_encap_t> fft;

        using index_t = tuple<size_t>;
        size_t translate_index(const index_t& index, const index_t& extends) const;

        typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType;
        Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate_matrix;
        Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict_matrix;
        Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> inject_matrix;



      public:

	//explicit SpectralTransfer(std::shared_ptr<fe_manager> FinEl, size_t nlevel);
	
        virtual void set_matrix(Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate, Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict, GFS& gfs1, GFS& gfs0);

	virtual void create(std::shared_ptr<fe_manager> FinEl, GFS& gfs1, GFS& gfs0);
	
        virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                      shared_ptr<typename TransferTraits::fine_encap_t> fine);

        virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);

      	virtual void restrict_u(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_t> coarse);
    };
  }  // ::pfasst::contrib
}  // ::pfasst

#include "spectral_transfer_impl.hpp"
//#include "pfasst/contrib/spectral_transfer_1d_impl.hpp"

#endif // _PFASST__TRANSFER__SPECTRAL_TRANSFER_1D_HPP_
