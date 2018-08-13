

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












const size_t GRID_LEVEL=1;


 class fe_manager{
	  
	  size_t n_elem;
	  size_t* n_dof;
	  size_t n_levels;
	  
          std::shared_ptr<GridType> grid;

	  std::shared_ptr<TransferOperatorAssembler<Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> >>> transfer;
	  std::shared_ptr<std::vector<MatrixType*>> transferMatrix;

		    
	  public:
	  
	   std::vector<std::shared_ptr<BasisFunction> > fe_basis; 
	   std::vector<std::shared_ptr<BasisFunction> > fe_basis_p; 
 
	    
	  fe_manager(const size_t nelements, size_t nlevels=1, MPI_Comm comm_x=MPI_COMM_WORLD, size_t base_order=1)
	  :fe_basis(nlevels), n_levels(nlevels)
	  {
	    n_elem=nelements;
	    n_dof = new size_t [nlevels];

	
	            const unsigned int dim = 2;
        Dune::FieldVector<typename GridType::ctype,dim> L;
        L[0]=1; L[1]=1;
        typename std::array<int,dim> s;
        std::fill(s.begin(), s.end(), nelements);
        std::bitset<dim> periodic;//(true, true);
        periodic[0]=true; periodic[1]=true;

        grid        = std::make_shared<GridType>(L,s,periodic,0);
    
	    

	    for (int i=0; i<nlevels; i++){
	      
	      grid->globalRefine(1);
              fe_basis[nlevels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i +1)); 
	      std::cout << i << " " << nlevels-i-1 << " *****  Ordnung " << fe_basis[nlevels-i-1]->size() << std::endl;              
	      n_dof[nlevels-i-1]    = fe_basis[nlevels-i-1]->size();

	    } 


	    std::cout << "***** 0 Ordnung " << fe_basis[0]->size() << std::endl;

	    std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
	    if(nlevels>1){ 
	      this->create_transfer();
	      
	    }

	    
	  }
	  
	  size_t get_ndofs(size_t i){return n_dof[i];}
	  size_t get_nelem(){return n_elem;}
	  std::shared_ptr<BasisFunction> get_basis(size_t i){return fe_basis[i];}
	  std::shared_ptr<GridType> get_grid(){return grid;}
	  std::shared_ptr<std::vector<MatrixType*>> get_transfer(){	   return transferMatrix;}
	  size_t get_nlevel() {return n_levels;}
	  
	  void create_transfer(){
	    transfer = std::make_shared<TransferOperatorAssembler<Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1>>>>(*grid);
	    transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	    for (int i=0; i< n_levels; i++){
	      transferMatrix->push_back(new MatrixType()); 
	    }
	    transfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	    
	    std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;

	  }
	  
	  
	  
};
