

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

//typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType;
//typedef Dune::BlockVector<Dune::FieldVector<double,1> > VectorType;

 class fe_manager{
	  
	  size_t n_elem;
	  size_t* n_dof;
	  size_t n_levels;

	  

	  //Dune::Functions::PQkNodalBasis<GridType::LeafGridView GridView,BASE_ORDER>;

          std::shared_ptr<GridType> grid;
	  
	  //std::shared_ptr<BasisFunction> basis;
	  //std::vector<BasisFunction> basis;
	  //std::shared_ptr<std::vector<std::shared_ptr<BasisFunction>>> basis; 
	  //std::vector<std::shared_ptr<BasisFunction> > fe_basis;

	  std::shared_ptr<TransferOperatorAssembler<Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> >>> transfer;
	  std::shared_ptr<std::vector<MatrixType*>> transferMatrix;
	  //MatrixType m1;
	  //std::vector<MatrixType> m;
		    
	  public:
	  
	   std::vector<std::shared_ptr<BasisFunction> > fe_basis; 
	   std::vector<std::shared_ptr<BasisFunction> > fe_basis_p; 
	  //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dune;
          //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > A_dune;  
	    
	  fe_manager(const size_t nelements, size_t nlevels=1, MPI_Comm comm_x=MPI_COMM_WORLD, size_t base_order=1)
	  :fe_basis(nlevels), n_levels(nlevels)
	  {
	    //Konstruktor
	    //hier wird das Gitter gebaut und die Basis des FE-Raums gewaehlt

	    n_elem=nelements;
	    n_dof = new size_t [nlevels];

	
	    const int DIMENSION=1;
	    //if(DIMENSION==2){
	      //Dune::FieldVector<double,DIMENSION> h = {1, 1};
	    //}
	    std::cout << DIMENSION << std::endl; 
	    
	    Dune::FieldVector<double,DIMENSION> h = {1};//-20, 20};
	    Dune::FieldVector<double,DIMENSION> hR = {200};
	    Dune::FieldVector<double,DIMENSION> hL = {-200};
	    array<int,DIMENSION> n;

	    std::fill(n.begin(), n.end(), nelements);

        /*int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

        	MPI_Comm comm_x, comm_t; 
	int myid, xcolor, tcolor;

	int space_num=2;
   	xcolor = (my_rank / space_num);
   	tcolor = my_rank % space_num;

   	MPI_Comm_split( MPI_COMM_WORLD, xcolor, my_rank, &comm_x );
   	MPI_Comm_split( MPI_COMM_WORLD, tcolor, my_rank, &comm_t );*/


#if HAVE_MPI
 	//grid = std::make_shared<GridType>(hL, hR, n, std::bitset<DIMENSION>{0ULL}, 1, comm_x); //overlap
 	this->grid = std::make_shared<GridType>(hL, hR, n, std::bitset<DIMENSION>{0ULL}, 1, comm_x); 
#else
        this->grid = std::make_shared<GridType>(hL, hR, n);
        //grid_global = std::make_shared<GridType>(hL, hR, n);
#endif	    
	    
	    //this->grid = std::make_shared<GridType>(hL, hR, n);
	    //this->grid = std::make_shared<GridType>(h, n);
	    //grid.globalRefine(0);
	    
	    
	    
	    //Dune::FieldVector<double,DIMENSION> h = {1};
	    
	      
	    //array<int,DIMENSION> n;
	    //std::fill(n.begin(), n.end(), nelements);

	    //this->grid  = std::make_shared<GridType>(h,n);
	    //this->basis = std::make_shared<std::vector<BasisFunction*>>(nlevels);
	    //std::vector<std::shared_ptr<BasisFunction>> test  
	    
	    
	    //std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
	    //std::cout << "***** Ordnung der Basis " << BASE_ORDER << std::endl;
	    
	    //this->basis = new std::vector<BasisFunction>(nlevels);
	    
	    for (int i=0; i<nlevels; i++){
	      
	      grid->globalRefine(1);//mlsdc 1
	      //GridType::LeafGridView gridView = grid->leafGridView();
	      //BasisFunction* b = new BasisFunction(gridView);
	      //std::cout << "***** groesse" << b->size() << std::endl;
	    
	      //basis->push_back(std::make_shared<BasisFunction>(gridView));	
	    
	      //std::cout << "***** groesse" << (&((&basis)[0]))->size() << std::endl;
	      //(*basis)[i]->gridView();
	   
	    
	      //while(basis.size() < nlevels)
              fe_basis[nlevels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i +1)); //grid->levelGridView(i));//gridView);
	      std::cout << i << " " << nlevels-i-1 << " *****  Ordnung " << fe_basis[nlevels-i-1]->size() << std::endl;              
	      n_dof[nlevels-i-1]    = fe_basis[nlevels-i-1]->size();

	    } 
	    //std::shared_ptr<BasisFunction> ruth = (*basis)[0];
	    //ruth->size();
	    //(*basis)[0]->gridView();

	    std::cout << "***** 0 Ordnung " << fe_basis[0]->size() << std::endl;
	    //std::cout << "***** 1 Ordnung " << fe_basis[1]->size() << std::endl;
	    //std::cout << "***** Ordnung " << fe_basis[1]->size() << std::endl;
	    //std::cout << "***** Ordnung " << fe_basis[2]->size() << std::endl;

	    std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
	    if(nlevels>1){ 
	      this->create_transfer();
	      //m.resize(nlevels);
	      
	    }
	    std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
	    //this->basis = std::make_shared<BasisFunction>(gridView);

	    //std::cout << "***** Basis erstellt mit " <<  fe_basis[0]->size() << " Elementen " << std::endl;

	    //n_dof=((fe_basis)[0])->size();
	    //std::cout << "***** Ordnung der Basis " << BASE_ORDER << std::endl;
	    //this->encap_factory()->set_size(basis->size());
	    //assembleProblem(((basis)[0]), A_dune, M_dune);
	    //std::cout << "***** Ordnung der Basis " << BASE_ORDER << std::endl;
	    
	    
	  }
	  
	  size_t get_ndofs(size_t i){return n_dof[i];}
	  size_t get_nelem(){return n_elem;}
	  std::shared_ptr<BasisFunction> get_basis(size_t i){return fe_basis[i];}
	  std::shared_ptr<GridType> get_grid(){return grid;}
	  //MatrixType get_transfer(size_t l){	    std::cout <<  "transfer rueckgabe" <<  std::endl; return *transferMatrix->at(0);}
	  std::shared_ptr<std::vector<MatrixType*>> get_transfer(){	   return transferMatrix;}
	  size_t get_nlevel() {return n_levels;}
	  
	  void create_transfer(){
	    std::cout << "im create transfer  neu " << std::endl;
	    transfer = std::make_shared<TransferOperatorAssembler<Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1>>>>(*grid);
	    transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	    	    std::cout << "im create transfer 1 " << std::endl;
	    for (int i=0; i< n_levels; i++){
	      transferMatrix->push_back(new MatrixType()); // hier nur referenz die evtl geloescht wird??
	    }
	    transfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	    
	    std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;
	    	    	    std::cout << "im create transfer vecvec " <<  (*vecvec->at(0)).M() << (*vecvec->at(0)).N() << (*vecvec->at(1)).M() << (*vecvec->at(1)).N() << std::endl;
	    std::cout <<  "transfer erzeugt groesse " << (*vecvec->at(1)).M() <<  std::endl;
	    for (int i=0; i< vecvec->at(0)->N(); i++){
	      for (int j=0; j< (*vecvec->at(0)).M(); j++){
		if(vecvec->at(0)->exists(i,j)){
		  //std::cout << ((*vecvec->at(0))[i][j]) << std::endl;
		}
	      }

            }
            //std::exit(0);
	  }
	  
	  
	  
};
