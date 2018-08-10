
//#include <dune/grid/yaspgrid.hh>
//#include "assemble.hpp"
#include <dune/fufem/assemblers/transferoperatorassembler_t.hh>
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









typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType2;
typedef Dune::BlockVector<Dune::FieldVector<double,1> > VectorType;

 class fe_manager{
	  
	  size_t n_elem;
	  size_t* n_dof;
	  size_t n_levels;

	  typedef Dune::YaspGrid<1> GridType; 
          typedef GridType::LevelGridView GridView;
	  using BasisFunction = Dune::Functions::PQkNodalBasis<GridView,1>;

          std::shared_ptr<GridType> grid;
	 
	  std::shared_ptr<TransferOperatorAssembler<Dune::YaspGrid<1>>> transfer;
	  std::shared_ptr<TransferOperatorAssembler_t<Dune::YaspGrid<1>>> transfer_t;
	  std::shared_ptr<std::vector<MatrixType2*>> transferMatrix;
	  std::shared_ptr<std::vector<MatrixType2*>> transferMatrix_t;

		    
	  public:
	  
	   std::vector<std::shared_ptr<BasisFunction> > fe_basis; 
	   std::vector<std::shared_ptr<BasisFunction> > fe_basis_p; 
	    
	  fe_manager(const size_t nelements, size_t nlevels=1, size_t base_order=1)
	  :fe_basis(nlevels), n_levels(nlevels)
	  {


	    n_elem=nelements;
	    n_dof = new size_t [nlevels];

	
	    const int DIMENSION=1;
	    Dune::FieldVector<double,DIMENSION> h = {1};
	    
	      
	    array<int,DIMENSION> n;
	    std::fill(n.begin(), n.end(), nelements);



#if HAVE_MPI
 	    this->grid = std::make_shared<GridType>(h,n, std::bitset<DIMENSION>{0ULL}, 1); //, MPI_COMM_SELF
#else
            this->grid = std::make_shared<GridType>(h,n);
#endif



	    
	    for (int i=0; i<nlevels; i++){
	      
	      grid->globalRefine((bool) i);

	      auto view = grid->levelGridView(i);
	      fe_basis[nlevels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); 
	      n_dof[nlevels-i-1]    = fe_basis[nlevels-i-1]->size();

	    } 


	    if(nlevels>1){ 
	      this->create_transfer();
	      
	    }
	    
	    
	  }
	  
	  size_t get_ndofs(size_t i){return n_dof[i];}
	  size_t get_nelem(){return n_elem;}
	  std::shared_ptr<BasisFunction> get_basis(size_t i){return fe_basis[i];}
	  std::shared_ptr<GridType> get_grid(){return grid;}
	  std::shared_ptr<std::vector<MatrixType2*>> get_transfer(){	   return transferMatrix;}
	  std::shared_ptr<std::vector<MatrixType2*>> get_transfer_t(){	   return transferMatrix_t;}
	  size_t get_nlevel() {return n_levels;}
	  
	  void create_transfer(){
	    transfer = std::make_shared<TransferOperatorAssembler<Dune::YaspGrid<1>>>(*grid);
	    transfer_t = std::make_shared<TransferOperatorAssembler_t<Dune::YaspGrid<1>>>(*grid);
	    transferMatrix = std::make_shared<std::vector<MatrixType2*>>();
	    transferMatrix_t = std::make_shared<std::vector<MatrixType2*>>();
	    for (int i=0; i< n_levels-1; i++){
	      transferMatrix->push_back(new MatrixType2());
 	      transferMatrix_t->push_back(new MatrixType2());
	    }
	    transfer->assembleMatrixHierarchy<MatrixType2>(*transferMatrix);
	    transfer_t->assembleMatrixHierarchy<MatrixType2>(*transferMatrix_t);
		/*		int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
		if(my_rank==0)
	    	for(int i=0; i<((*transferMatrix)[0][0]).N(); i++){
		for(int j=0; j<((*transferMatrix)[0][0]).M(); j++){ 
			if ((*transferMatrix)[0][0].exists(i,j)) {
				std::cout << (*transferMatrix)[0][0][i][j][0][0]<< " ";
			}else { 
				std::cout << 0;
			} 
		}std::cout << 0 << std::endl;
		
	}std::cout << std::endl;std::exit(0);*/

	  }
	  
	  
	  
};
