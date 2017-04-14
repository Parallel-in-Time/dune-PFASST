
//#include <dune/grid/yaspgrid.hh>
//#include "assemble.hpp"
#include <dune/fufem/assemblers/transferoperatorassembler.hh>


#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>


#include <dune/functions/functionspacebases/interpolate.hh>

#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>


#include <dune/fufem/assemblers/basisinterpolationmatrixassembler.hh>










typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType;
typedef Dune::BlockVector<Dune::FieldVector<double,1> > VectorType;

 class fe_manager{
	  
	  size_t n_elem;
	  size_t* n_dof;
	  size_t n_levels;

	  typedef Dune::YaspGrid<1> GridType; 
      typedef GridType::LeafGridView GridView;
	  //using BasisFunction = Dune::Functions::PQkNodalBasis<GridView,1>;
      
      typedef DuneFunctionsBasis<Dune::Functions::PQkNodalBasis<GridType::LeafGridView,1>> B1;
      typedef DuneFunctionsBasis<Dune::Functions::PQkNodalBasis<GridType::LeafGridView,2>> B2;

      

      std::shared_ptr<GridType> grid;
	  


	  std::shared_ptr<TransferOperatorAssembler<Dune::YaspGrid<1>>> transfer;
	  std::shared_ptr<std::vector<MatrixType*>> transferMatrix;

		    
	  public:
	  
	   std::shared_ptr<B1>  fe_basis1; 
	   std::shared_ptr<B2>  fe_basis2; 

	    
	  fe_manager(const size_t nelements, size_t nlevels=1, size_t base_order=1)
	  : n_levels(nlevels)
	  {


	    n_elem=nelements;
	    n_dof = new size_t [nlevels];

	
	    const int DIMENSION=1;

	    Dune::FieldVector<double,DIMENSION> h = {1};
	    
	      
	    array<int,DIMENSION> n;
	    std::fill(n.begin(), n.end(), nelements);

	    this->grid  = std::make_shared<GridType>(h,n);

        GridType::LeafGridView gridView = grid->leafGridView();
	    
        fe_basis2 = std::make_shared<B2>(gridView); 
	    n_dof[0]    = fe_basis2->size();
        
	    fe_basis1 = std::make_shared<B1>(gridView); 
	    n_dof[1]    = fe_basis1->size();
        

	    
	    
	    if(nlevels>1){ 
	      //this->create_transfer();
	      
	    }

	    
	    
	  }
	  
	  size_t get_ndofs(size_t i){return n_dof[i];}
	  size_t get_nelem(){return n_elem;}
	  std::shared_ptr<B2> get_basis0(){
          return fe_basis2;        
      }
      std::shared_ptr<B1> get_basis1(){
          return fe_basis1;          
      }
    
	  std::shared_ptr<GridType> get_grid(){return grid;}

	  std::shared_ptr<std::vector<MatrixType*>> get_transfer(){	   return transferMatrix;}
	  size_t get_nlevel() {return n_levels;}
	  
	  /*void create_transfer(){
	    transfer = std::make_shared<TransferOperatorAssembler<Dune::YaspGrid<1>>>(*grid);
	    transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	    for (int i=0; i< n_levels-1; i++){
	      transferMatrix->push_back(new MatrixType()); 
	    }
	    transfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	    
	    std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;

	    for (int i=0; i< vecvec->at(0)->N(); i++){
	      for (int j=0; j< (*vecvec->at(0)).M(); j++){
		if(vecvec->at(0)->exists(i,j)){

		}
	      }

        }
	  }*/
	  
	  
	  
};
