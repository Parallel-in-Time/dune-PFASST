#include <mpi.h>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/densematrix.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/typetree/utility.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>

#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/istlbackend.hh>
#include <dune/fufem/formatstring.hh>
#include <dune/fufem/assemblers/transferoperatorassembler.hh>

const int BASE_ORDER=1;
const int DIMENSION=1;
const int NR_COMP=1;
const int nelements=128;
const int n_levels=2;
const double _nu=1;
const double _n=1;
const double dt=0.1;

typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,NR_COMP,NR_COMP> > MatrixType;
typedef Dune::BlockVector<Dune::FieldVector<double,NR_COMP> > VectorType;
typedef Dune::FieldVector<double,NR_COMP>   InVectorType;
typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
typedef GridType::LevelGridView GridView;

using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
using namespace std;

//calculates Mu + dt(_nu (wu^(_n+1) -Mu) -Au)-rhs
void evaluate_f(VectorType &f,
	  const VectorType u,
          const VectorType rhs,
	  const MatrixType M_dune,
          const MatrixType A_dune,
	  const VectorType w){

        f*=0;
        for (int i=0; i<u.size(); ++i) f[i]= pow(u[i], _n+1) * (w)[i];	
        M_dune.mmv(u, f);
        f *= (_nu*_nu);
        A_dune.mmv(u,f);
        f *= dt;
        M_dune.umv(u,f);
        f -=rhs;
}

//calculate derivative of evaluate_f						
void evaluate_df(MatrixType &df,
	   const VectorType u,
	   const MatrixType M_dune,
           const MatrixType A_dune, 
           const VectorType w){

        for (int i=0; i<df.N(); ++i) df[i][i]= (_nu*_nu)*(_n+1) * pow(u[i], _n) * ((double) (w)[i]);	
        df.axpy((-_nu*_nu), M_dune);
        df-=A_dune;
        df*=dt;
        df+=M_dune;
} 

int main(int argc, char** argv) {

  	MPI_Init(&argc, &argv);
        //typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
        //typedef GridType::LevelGridView GridView;
	//using BasisFunction = Dune::Functions::PQkNodalBasis<GridView,1>;


        std::shared_ptr<GridType> grid;   
	std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); //Basen auf den verschiedenen Gitterebenen 


 	Dune::FieldVector<double,DIMENSION> hR = {200};
 	Dune::FieldVector<double,DIMENSION> hL = {-200};
 	array<int,DIMENSION> n;
 
 	std::fill(n.begin(), n.end(), nelements);
  	    
#if HAVE_MPI
 	grid = std::make_shared<GridType>(hL, hR, n, std::bitset<DIMENSION>{0ULL}, 1, MPI_COMM_SELF);
#else
        grid = std::make_shared<GridType>(hL, hR, n);
#endif
	    
	for (int i=0; i<n_levels; i++){	      
	      grid->globalRefine((bool) i);
	      auto view = grid->levelGridView(i);
              fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); 
	}
	BasisFunction basis = *fe_basis[0];//n_levels-1];
	VectorType coarse(fe_basis[n_levels-1]->size()); 




	//std::exit(0);
      	MatrixType M_dune;
      	MatrixType A_dune;

        auto sBackend = Dune::Fufem::istlMatrixBackend(A_dune);
        auto mBackend = Dune::Fufem::istlMatrixBackend(M_dune);

        using Assembler = Dune::Fufem::DuneFunctionsOperatorAssembler<BasisFunction, BasisFunction>;
        auto assembler = Assembler{basis, basis};

        using FiniteElement = std::decay_t<decltype(basis.localView().tree().finiteElement())>;

        auto vintageLaplace = LaplaceAssembler<GridType,FiniteElement, FiniteElement>();
        auto vintageMass = MassAssembler<GridType,FiniteElement, FiniteElement>();

        auto localAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageLaplace.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
        };
        auto localMassAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageMass.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
        };

        assembler.assembleBulk(sBackend, localAssembler);
        assembler.assembleBulk(mBackend, localMassAssembler);
                
        A_dune*=-1;
    
        std::shared_ptr<VectorType> w; 
        w = std::make_shared<VectorType>(M_dune.M());
        
        for(int j=0; j<M_dune.M(); j++){(*w)[j]=0;}

        for(int i=0; i<M_dune.M(); i++){
        	for(int j=0; j<M_dune.M(); j++){
                        if(M_dune.exists(i,j)) (*w)[i][0]= ((double) (*w)[i][0]) + ((double) M_dune[i][j][0][0]);
        	}
	}

	VectorType residuum(M_dune.M()), u(M_dune.M()), rhs(M_dune.M()), f(M_dune.M()), exact(M_dune.M()), newton_rhs(M_dune.M());


    
    

	//const auto dim = 1; //SweeperTrait::DIM;
	double t=0;
        double l0 = _nu;
        double l1 = l0/2. *(pow((1+_n/2.), 1/2.) + pow((1+ _n/2.), -1/2.) );
        double d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
        auto initial = [l0, l1, _n, d, t](const Dune::FieldVector<double,DIMENSION>&x){ 
                return pow((1 + (pow(2, _n/2.)-1 )* exp(-(_n/2.)*d*(x+2*l1*t)) ), -2./_n);
        };  
        interpolate(basis, u, initial);
        
        double t_end=0.1;
        auto exact_solution = [l0, l1, _n, d, t_end](const Dune::FieldVector<double,DIMENSION> &x){ 
                return pow((1 + (pow(2, _n/2.)-1 )* exp(-(_n/2.)*d*(x+2*l1*t_end)) ), -2./_n);
        };  
        interpolate(basis, exact, exact_solution);
        M_dune.mv(u, rhs);
        
    	for(int i=0; i<rhs.size(); i++){std::cout << "***** Anfangswert " << u[i] << std::endl;}
    
	
	for (int i=0; i< 5 ;i++){
	    
            MatrixType df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(M_dune); 
            evaluate_f(f, u, rhs,M_dune, A_dune, *w);
            evaluate_df(df, u, M_dune, A_dune, *w);
            df.mv(u, newton_rhs);
            newton_rhs -= f;
          
            Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df);
	  
	    df[0][0]=1;
            df[0][1]=0;

            df[df.N()-1][df.M()-1]=1;
            df[df.N()-1][df.M()-2]=0;

            Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df,1.0);	  
            Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-16, 
                              5000,    
                              0);    
          
            Dune::InverseOperatorResult statistics ;
            newton_rhs[0]=0; newton_rhs[rhs.size()-1]=1;
            cg.apply(u, newton_rhs , statistics ); //rhs is not const

            evaluate_f(f, u, rhs, M_dune, A_dune, *w);          
            std::cout << i << " residuumsnorm von f(u) " << f.infinity_norm() << std::endl;  
            if(f.infinity_norm()<1e-10){   std::cout << "genauigkeit erreicht " << i << std::endl;      break;} 
	}


	for(int i=0; i< M_dune.N(); i++){
		for(int j=0; j< M_dune.M(); j++){
		if(M_dune.exists(i,j)){ std::cout << M_dune[i][j]<< " ";}else{std::cout <<" ";}
		}std::cout << std::endl;
	}

	for(int i=0; i<u.size(); i++){std::cout << "***** Loesung " << u[i] << " exact "<< exact[i] << std::endl;}

	std::shared_ptr<TransferOperatorAssembler<GridType>> transfer;
	std::shared_ptr<std::vector<MatrixType*>> transferMatrix;
	  
	transfer = std::make_shared<TransferOperatorAssembler<GridType>>(*grid);
	transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	for (int i=0; i< n_levels-1; i++){
	      transferMatrix->push_back(new MatrixType()); 
	}
	transfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	MPI_Barrier(MPI_COMM_WORLD);
	
	std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;


	
	MatrixType interpolate_matrix= *(vecvec->at(0));
	MatrixType restrict_matrix   = *(vecvec->at(0));

	//std::cout << interpolate_matrix.N() <<" "<< interpolate_matrix.M() <<" "<< restrict_matrix.M() <<" "<< restrict_matrix.N() << std::endl;
	//std::cout << fe_basis[n_levels-2]->size() << " " << fe_basis[n_levels-1]->size() << std::endl;
	//std::exit(0);
	for (int i=0; i< restrict_matrix.N(); i++){
		for (int j=0; j< restrict_matrix.M(); j++){
			if(restrict_matrix.exists(i,j)){	
		  		if (restrict_matrix[i][j]==0.5 ) restrict_matrix[i][j]=0;
			}

	      	}
	}
	
	for (int i=0; i< interpolate_matrix.N(); i++){
	      for (int j=0; j< interpolate_matrix.M(); j++){
		if(interpolate_matrix.exists(i,j)){
		  std::cout << (interpolate_matrix[i][j]) << std::endl;
		}
	      }

        }
	
	for (int i=0; i< restrict_matrix.N(); i++){
	      for (int j=0; j< restrict_matrix.M(); j++){
		if(restrict_matrix.exists(i,j)){
		  std::cout << (restrict_matrix[i][j]) << std::endl;
		}
	      }

        }

	restrict_matrix.mtv(u, coarse);
	interpolate_matrix.mtv(u, coarse);
	
	for (int i=0; i<coarse.size(); i++)  std::cout <<i << " " << coarse[i] << " test" << std::endl;	
	interpolate_matrix.mv(coarse, u);

	
	MPI_Finalize();
}




