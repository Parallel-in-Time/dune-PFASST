
#include <config.h>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <cstring>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/densematrix.hh>
#include<dune/common/parallel/mpihelper.hh>

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
#include <dune/istl/io.hh>
#include<dune/istl/matrixmarket.hh>
#include<dune/istl/matrixredistribute.hh>
#include <dune/istl/schwarz.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/typetree/utility.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#if USE_DG
#  include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#else
#  include <dune/functions/functionspacebases/pqknodalbasis.hh>
#endif
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/istlbackend.hh>
#include <dune/fufem/formatstring.hh>
#include <dune/fufem/assemblers/transferoperatorassembler.hh>

#include<dune/istl/paamg/pinfo.hh>
#include<dune/istl/paamg/graph.hh>
#if USE_DG
#  include <dune/parmg/test/dglaplacematrix.hh>
#else
#  include <dune/parmg/test/laplacematrix.hh>
#endif
#include <dune/parmg/iterationstep/lambdastep.hh>
#include <dune/parmg/iterationstep/multigrid.hh>
#include <dune/parmg/iterationstep/multigridstep.hh>
#include <dune/parmg/norms/normadapter.hh>
#include <dune/parmg/parallel/communicationp1.hh>
#include <dune/parmg/parallel/communicationdg.hh>
#include <dune/parmg/parallel/datahandle.hh>
#include <dune/parmg/parallel/dofmap.hh>
#include <dune/parmg/parallel/globaldofindex.hh>
#include <dune/parmg/parallel/istlcommunication.hh>
#include <dune/parmg/parallel/matrixalgebra.hh>
#include <dune/parmg/parallel/vectoralgebra.hh>
#include <dune/parmg/parallel/redistributematrix.hh>
#include <dune/parmg/parallel/redistributevector.hh>
#include <dune/parmg/parallel/parallelenergyfunctional.hh>
#include <dune/parmg/parallel/parallelenergynorm.hh>
#include <dune/parmg/parallel/restrictmatrix.hh>
#include <dune/parmg/solvers/coarsesuperlusolver.hh>
#include <dune/parmg/solvers/linesearch.hh>
#include <dune/parmg/solvers/directionsearch.hh>
#include <dune/parmg/iterationstep/multigridsetup.hh>
#include <dune/parmg/iterationstep/parallelprojectedgs.hh>

#include <dune/solvers/iterationsteps/blockgssteps.hh>
#include <dune/solvers/norms/energynorm.hh>
#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/transferoperators/compressedmultigridtransfer.hh>




const int BASE_ORDER=1;
const int DIMENSION=1;
const int NR_COMP=1;
const int nelements=128;
const int n_levels=3;
const double _nu=1;
const double _n=1;
const double dt=0.1;

typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,NR_COMP,NR_COMP> > MatrixType;
typedef Dune::BlockVector<Dune::FieldVector<double,NR_COMP> > VectorType;
typedef Dune::FieldVector<double,NR_COMP>   InVectorType;
typedef Dune::YaspGrid<DIMENSION,Dune::EquidistantOffsetCoordinates<double, DIMENSION> > GridType; 
typedef GridType::LevelGridView GridView;


using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
using FEBasis = BasisFunction; 
using namespace std;
using namespace Dune;

using namespace Dune;
using namespace Dune::ParMG;


//calculates Mu + dt(_nu (wu^(_n+1) -Mu) -Au)-rhs
void evaluate_f(VectorType &f,
	  const VectorType& u,
          const VectorType& rhs,
	  const MatrixType& M_dune,
          const MatrixType& A_dune,
	  const VectorType& w){

        f*=0;
        for (int i=0; i<u.size(); ++i) f[i]= pow(u[i], _n+1) * (w)[i];	
        M_dune.mmv(u, f);
        f *= (_nu*_nu);
        A_dune.umv(u,f);
        f *= dt;
        M_dune.umv(u,f);
        f -=rhs;
}

//calculate derivative of evaluate_f						
void evaluate_df(MatrixType &df,
	   const VectorType& u,
	   const MatrixType& M_dune,
           const MatrixType& A_dune, 
           const VectorType& w){

        df*=0;
        for (int i=0; i<df.N(); ++i) df[i][i]= (_nu*_nu)*(_n+1) * pow(u[i], _n) * ((double) (w)[i]);	
        df.axpy((-_nu*_nu), M_dune);
        df+=A_dune;
        df*=dt;
        df+=M_dune;
} 

int main(int argc, char** argv) {

  	MPI_Init(&argc, &argv);
	
	int rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        

        std::shared_ptr<GridType> grid;   

	std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); //Basen auf den verschiedenen Gitterebenen 

 	Dune::FieldVector<double,DIMENSION> hR = {200};
 	Dune::FieldVector<double,DIMENSION> hL = {-200};

 	array<int,DIMENSION> n;

 	std::fill(n.begin(), n.end(), nelements);

  	    
#if HAVE_MPI
 	grid = std::make_shared<GridType>(hL, hR, n, std::bitset<DIMENSION>{0ULL}, 1, MPI_COMM_WORLD); //overlap
#else
        grid = std::make_shared<GridType>(hL, hR, n);
#endif
	
  
	//grid->globalRefine(1);	 
	for (int i=0; i<n_levels; i++){	      
	      grid->globalRefine(1);
              fe_basis[i] = std::make_shared<BasisFunction>(grid->levelGridView(i+1)); 
	}

	int coarser_level=2;
	BasisFunction basis = *fe_basis[coarser_level];

	      ;
	//VectorType coarse(fe_basis[n_levels-1]->size()); 




	
	//setup stiffness and massmatrix use lumping for nonlinear part
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
    
        std::shared_ptr<VectorType> w; 
        w = std::make_shared<VectorType>(M_dune.M());

        
        for(int j=0; j<M_dune.M(); j++){(*w)[j]=0;}

        for(int i=0; i<M_dune.M(); i++){
        	for(int j=0; j<M_dune.M(); j++){
                        if(M_dune.exists(i,j)) (*w)[i][0]= ((double) (*w)[i][0]) + ((double) M_dune[i][j][0][0]);
        	}
	}



	
	VectorType residuum(M_dune.M()), u(M_dune.M()), rhs(M_dune.M()), f(M_dune.M()), exact(M_dune.M()), newton_rhs(M_dune.M());

    
    

	//set initial condition 
	double t=0;
        double l0 = _nu;
        double l1 = l0/2. *(pow((1+_n/2.), 1/2.) + pow((1+ _n/2.), -1/2.) );
        double d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
        auto initial = [l0, l1, _n, d, t](const Dune::FieldVector<double,1>&x){ 
                return pow((1 + (pow(2, _n/2.)-1 )* exp(-(_n/2.)*d*(x+2*l1*t)) ), -2./_n);
        };  
        interpolate(basis, u, initial);

        double t_end=0.1;
        auto exact_solution = [l0, l1, _n, d, t_end](const Dune::FieldVector<double,1> &x){ 
                return pow((1 + (pow(2, _n/2.)-1 )* exp(-(_n/2.)*d*(x+2*l1*t_end)) ), -2./_n);
        };  
        interpolate(basis, exact, exact_solution);

        M_dune.mv(u, rhs);


	

	//start newton to solve the resulting nonlinear system
	for (int i=0; i< 10 ;i++){
	    
	    //setup jacobi-matrix df and newton_rhs
            MatrixType df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(M_dune); 

            evaluate_f(f, u, rhs,M_dune, A_dune, *w);

            evaluate_df(df, u, M_dune, A_dune, *w);


            
            df.mv(u, newton_rhs);
            newton_rhs -= f;






///////////////////////////////////////////////////////////////////////////////////////

	//here we use parmg solver to solve the local system (df * u = newton_rhs) for the unknown u
  	auto &x=u;
  	using MGSetup = Dune::ParMG::ParallelMultiGridSetup< FEBasis, MatrixType, VectorType >;
  	
  	std::cout << "vor erster grid aktion" << std::endl;
  	MPI_Barrier(MPI_COMM_WORLD);
  	MGSetup mgSetup{*grid,coarser_level+1}; //0 grobgitter
  	std::cout << "nach grid setup " << std::endl;
  	auto gridView = mgSetup.bases_.back().gridView();
 	using MG = Dune::ParMG::Multigrid<VectorType>;
  	MG mg;
  
    	using namespace Dune::ParMG;
    	auto& levelOp = mgSetup.levelOps_;
    	auto df_pointer = std::make_shared<MatrixType>(df);	
    	mgSetup.matrix(df_pointer);
    	auto fineIgnore = std::make_shared< BitSetVector<1> >(u.size());
    	for (std::size_t i = 0; i < u.size(); ++i){
      		(*fineIgnore)[i] = false;
      		if(rank==num_pro-1 && i== u.size()-1) (*fineIgnore)[i] = true;
      	}
    	mgSetup.ignore(fineIgnore);
    	mgSetup.setupLevelOps();
    	double dampening =1.0;
    	mgSetup.setupSmoother(dampening);
    	bool enableCoarseCorrection=true;
    	if (enableCoarseCorrection)
      		mgSetup.setupCoarseSuperLUSolver();
    	else
      		mgSetup.setupCoarseNullSolver();
    	mg.levelOperations(levelOp);
    	mg.coarseSolver(mgSetup.coarseSolver());
    	levelOp.back().maybeRestrictToMaster(newton_rhs);
    	std::function<void(VectorType&)> collect = Dune::ParMG::makeCollect<VectorType>(*mgSetup.comms_.back());
    	std::function<void(VectorType&)> restrictToMaster = [op=levelOp.back()](VectorType& x) { op.maybeRestrictToMaster(x); };
    	auto energyFunctional = Dune::ParMG::makeParallelEnergyFunctional(
      		*df_pointer,
      		newton_rhs,
      		gridView.grid().comm(),
      		//collect
      		restrictToMaster
      	);
    	auto energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(*df_pointer, restrictToMaster, gridView.grid().comm());
    	levelOp.back().maybeCopyFromMaster(x);
    	double tol = 1e-13;
    	auto realIterationStep = [&](VectorType& x) {
      		auto b = newton_rhs;
      		mg.apply(x, b);
    	};
  	auto& feBasis = mgSetup.bases_.back();
    	auto solverNorm = std::make_shared< NormAdapter<VectorType> >(energyNorm);
	auto iterationStep = std::make_shared< LambdaStep<VectorType> >(realIterationStep, x);
	int steps = 100;
    	auto solver = Dune::Solvers::LoopSolver<VectorType>(iterationStep, steps, tol, solverNorm, NumProc::FULL);
	solver.preprocess();
    	solver.solve();
    

//////////////////////////////////////////////////////////////////////////////////////






	    for (int k=0; k<u.size(); k++)  if (rank==0) std::cout << "rank 0 Newton " << i << " Index " << k<< " local u " << u[k] << " " << std::endl;
	    MPI_Barrier(MPI_COMM_WORLD);
	    for (int k=0; k<u.size(); k++)  if (rank==num_pro-1) std::cout << "rank 2 Newton " << i << " Index " << k<< " local u " << u[k] << " " << std::endl;
	    MPI_Barrier(MPI_COMM_WORLD);


            evaluate_f(f, u, rhs, M_dune, A_dune, *w); 
                     
                     
            using Norm =  EnergyNorm<MatrixType,VectorType>;
	    auto norm = Norm(A_dune);   
	    
	    
	    auto parallel_energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(A_dune, restrictToMaster, gridView.grid().comm());
	    
 
            std::cout << i << " residuumsnorm von f(u) " << parallel_energyNorm(f) << std::endl;  
            
            //if(rank==2) for(int i=0; i<rhs.size(); i++){std::cout << "0 *****  " << u[i] << std::endl;}	
            //if(f_global.infinity_norm()<1e-10){   std::cout << "genauigkeit erreicht global " << i << std::endl;      break;} 	  
            //if (rank==1) if(f.infinity_norm()<1e-10){   std::cout << "genauigkeit erreicht local" << i << std::endl;      break;} 	  

	    

	}
	
	
	



	MPI_Finalize();
}




