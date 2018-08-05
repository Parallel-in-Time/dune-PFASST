
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
const int nelements=16000;
const int n_levels=2;
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

//solve frame propagation problem
int main(int argc, char** argv) {

  	MPI_Init(&argc, &argv);
	
	int rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        
	MPI_Comm comm_x, comm_t; 
	int myid, xcolor, tcolor;

	int space_num=2;
   	xcolor = (rank / space_num);
   	tcolor = rank % space_num;

   	MPI_Comm_split( MPI_COMM_WORLD, xcolor, rank, &comm_x );
   	MPI_Comm_split( MPI_COMM_WORLD, tcolor, rank, &comm_t );


        //pfasst.communicator() = std::make_shared<CommunicatorType>(comm_t);
	
	
	
	
	//setup the grid and basis functions //sequential vs parallel (*_global) setup
        std::shared_ptr<GridType> grid;   
        std::shared_ptr<GridType> grid_global;   
	std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); //Basen auf den verschiedenen Gitterebenen 
	std::vector<std::shared_ptr<BasisFunction> > fe_basis_global(n_levels);  

 	Dune::FieldVector<double,DIMENSION> hR = {200};
 	Dune::FieldVector<double,DIMENSION> hL = {-200};
 	Dune::FieldVector<double,DIMENSION> hR_global = {200};
 	Dune::FieldVector<double,DIMENSION> hL_global = {-200};
 	array<int,DIMENSION> n;
 	array<int,DIMENSION> n_global;
 	
 	std::fill(n.begin(), n.end(), nelements);
  	std::fill(n_global.begin(), n_global.end(), nelements);    
  	    
  	//MPI_Comm comm = MPI_COMM_WORLD;    
  	
  	/*Dune::FieldVector<typename GridType::ctype,DIMENSION> L;
        L[0]=1; L[1]=1;
        typename std::array<int,DIMENSION> s;
        std::fill(s.begin(), s.end(), nelements);
        std::bitset<DIMENSION> periodic;//(true, true);
        periodic[0]=false; periodic[1]=false;*/


  	
  	
  	//auto gridptr = Dune::StructuredGridFactory<GridType>::createCubeGrid({0,0},{1,1},{{2,2}});
#if HAVE_MPI
        //grid = std::make_shared<GridType>(L,s,periodic,1, comm_x);
        //grid_global = std::make_shared<GridType>(L,s,periodic,0,MPI_COMM_SELF);
 	grid = std::make_shared<GridType>(hL, hR, n, std::bitset<DIMENSION>{0ULL}, 1, MPI_COMM_WORLD); //overlap
 	grid_global = std::make_shared<GridType>(hL_global, hR_global, n_global, std::bitset<DIMENSION>{0ULL}, 0, MPI_COMM_SELF); 
#else
        grid = std::make_shared<GridType>(hL, hR, n);
        grid_global = std::make_shared<GridType>(hL, hR, n);
        //grid = std::make_shared<GridType>(L,s,periodic);
        //grid_global = std::make_shared<GridType>(L,s,periodic);
#endif
	
	int level_number=0;  
	grid->globalRefine(level_number);
	grid_global->globalRefine(level_number);   
	 
	for (int i=0; i<n_levels; i++){	      
	      grid->globalRefine((bool) i);
	      grid_global->globalRefine((bool) i);
	      auto view = grid->levelGridView(i);
	      auto view_global = grid->levelGridView(i);
              fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); 
              fe_basis_global[n_levels-i-1] = std::make_shared<BasisFunction>(grid_global->levelGridView(i)); 
	}
	BasisFunction basis = *fe_basis[0];//n_levels-1];
	//if(rank==0)std::cout << fe_basis[0]->size() << std::endl;
	//if(rank==0)std::cout << fe_basis[1]->size() << std::endl;
	//if(rank==0)std::cout << fe_basis[n_levels-1]->size() << std::endl;	
	//std::exit(0);

	BasisFunction basis_global = *fe_basis_global[0];//n_levels-1];
	VectorType coarse(fe_basis[n_levels-1]->size()); 
	VectorType coarse_global(fe_basis_global[n_levels-1]->size()); 


	std::cout << "globale Basis hat " << basis_global.size() << " Elemente und lokale " << basis.size()  << std::endl;
	
	//setup stiffness and massmatrix use lumping for nonlinear part
      	MatrixType M_dune;
      	MatrixType A_dune;
      	MatrixType M_dune_global;
      	MatrixType A_dune_global;

        auto sBackend = Dune::Fufem::istlMatrixBackend(A_dune);
        auto mBackend = Dune::Fufem::istlMatrixBackend(M_dune);
        auto sBackend_global = Dune::Fufem::istlMatrixBackend(A_dune_global);
        auto mBackend_global = Dune::Fufem::istlMatrixBackend(M_dune_global);

        using Assembler = Dune::Fufem::DuneFunctionsOperatorAssembler<BasisFunction, BasisFunction>;
        auto assembler = Assembler{basis, basis};
        auto assembler_global = Assembler{basis_global, basis_global};

        using FiniteElement = std::decay_t<decltype(basis.localView().tree().finiteElement())>;
        using FiniteElement_global = std::decay_t<decltype(basis_global.localView().tree().finiteElement())>;

        auto vintageLaplace = LaplaceAssembler<GridType,FiniteElement, FiniteElement>();
        auto vintageMass = MassAssembler<GridType,FiniteElement, FiniteElement>();
        auto vintageLaplace_global = LaplaceAssembler<GridType,FiniteElement_global, FiniteElement_global>();
        auto vintageMass_global = MassAssembler<GridType,FiniteElement_global, FiniteElement_global>();

        auto localAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageLaplace.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
        };
        auto localMassAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageMass.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
        };

        auto localAssembler_global = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageLaplace_global.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
        };
        auto localMassAssembler_global = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
                    vintageMass_global.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
        };


        assembler_global.assembleBulk(sBackend_global, localAssembler);
        assembler_global.assembleBulk(mBackend_global, localMassAssembler);

        assembler.assembleBulk(sBackend, localAssembler);
        assembler.assembleBulk(mBackend, localMassAssembler);         
                
        //this must be changed to the system otherwise the norm does not make sence!!!                
        //A_dune*=-1;
        //A_dune_global *=-1;
        //
        

       	
       	/*for (int i=0; i< A_dune.N(); i++){
		for (int j=0; j< A_dune.M(); j++){
			if(A_dune.exists(i,j)){	
		  		if(rank==0) std::cout << A_dune[i][j] << " ";
			}else{if(rank==0) std::cout << "0 ";}

	}
	
	
	for (int i=0; i< A_dune_global.N(); i++){
		for (int j=0; j< A_dune_global.M(); j++){
			if(A_dune_global.exists(i,j)){	
		  		if(rank==0) std::cout << A_dune_global[i][j] << " ";
			}else{if(rank==0) std::cout << "0 ";}

	      	}if (rank ==0) std::cout << std::endl;
	}
	std::cout << "nach der ausgabe " << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);*/
    
        std::shared_ptr<VectorType> w; 
        w = std::make_shared<VectorType>(M_dune.M());
        
        std::shared_ptr<VectorType> w_global; 
        w_global = std::make_shared<VectorType>(M_dune_global.M());
        
        for(int j=0; j<M_dune.M(); j++){(*w)[j]=0;}

        for(int i=0; i<M_dune.M(); i++){
        	for(int j=0; j<M_dune.M(); j++){
                        if(M_dune.exists(i,j)) (*w)[i][0]= ((double) (*w)[i][0]) + ((double) M_dune[i][j][0][0]);
        	}
	}

        for(int j=0; j<M_dune_global.M(); j++){(*w_global)[j]=0;}

        for(int i=0; i<M_dune_global.M(); i++){
        	for(int j=0; j<M_dune_global.M(); j++){
                        if(M_dune_global.exists(i,j)) (*w_global)[i][0]= ((double) (*w_global)[i][0]) + ((double) M_dune_global[i][j][0][0]);
        	}
	}

	
	VectorType residuum(M_dune.M()), u(M_dune.M()), u_neu(M_dune.M()), rhs(M_dune.M()), f(M_dune.M()), exact(M_dune.M()), newton_rhs(M_dune.M());
	VectorType residuum_global(M_dune_global.M()), u_global(M_dune_global.M()),u_neu_global(M_dune_global.M()), rhs_global(M_dune_global.M()), f_global(M_dune_global.M()), exact_global(M_dune_global.M()), newton_rhs_global(M_dune_global.M());

    
    

	//set initial condition 
	double t=0;
        double l0 = _nu;
        double l1 = l0/2. *(pow((1+_n/2.), 1/2.) + pow((1+ _n/2.), -1/2.) );
        double d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
        auto initial = [l0, l1, _n, d, t](const Dune::FieldVector<double,1>&x){ 
                return pow((1 + (pow(2, _n/2.)-1 )* exp(-(_n/2.)*d*(x+2*l1*t)) ), -2./_n);
        };  
        interpolate(basis, u, initial);
        interpolate(basis, u_neu, initial);
        interpolate(basis_global, u_global, initial);
        interpolate(basis_global, u_neu_global, initial);
        double t_end=0.1;
        auto exact_solution = [l0, l1, _n, d, t_end](const Dune::FieldVector<double,1> &x){ 
                return pow((1 + (pow(2, _n/2.)-1 )* exp(-(_n/2.)*d*(x+2*l1*t_end)) ), -2./_n);
        };  
        interpolate(basis, exact, exact_solution);
        interpolate(basis_global, exact_global, exact_solution);
        
        
        //ScalarRandomDiscFunction<dim> initialPhaseField(29007122, 50, 0.1, 0.5, -1.0, 1.0);
        /*auto phaseField = [&](auto x) {Dune::FieldVector<double, 1> y=1; return y;};
        Dune::Functions::interpolate(basis, u, phaseField);
        Dune::Functions::interpolate(basis_global, u_global, phaseField);
        Dune::Functions::interpolate(basis, exact, phaseField);
        Dune::Functions::interpolate(basis_global, exact_global, phaseField);*/
        
        M_dune.mv(u, rhs);
        M_dune_global.mv(u_global, rhs_global);
        

        /*if(rank==0) for(int i=0; i<rhs_global.size(); i++){std::cout << "*****  " << rhs_global[i] << std::endl;}
    	if(rank==0) for(int i=0; i<rhs.size(); i++){std::cout << "0 *****  " << rhs[i] << std::endl;}	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==1) for(int i=0; i<rhs.size(); i++){std::cout << "1 *****  " << rhs[i] << std::endl;}
	
	MPI_Barrier(MPI_COMM_WORLD);
    	//std::exit(0);*/
	

	//start newton to solve the resulting nonlinear system
	//for(int j=0; j < 100; j++)
	
	for (int i=0; i< 1000 ;i++){
	    
	    //setup jacobi-matrix df and newton_rhs
            MatrixType df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(M_dune); 
            MatrixType df_global = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(M_dune_global); 
            evaluate_f(f, u, rhs,M_dune, A_dune, *w);
            evaluate_f(f_global, u_global, rhs_global,M_dune_global, A_dune_global, *w_global);
            evaluate_df(df, u, M_dune, A_dune, *w);
            evaluate_df(df_global, u_global, M_dune_global, A_dune_global, *w_global);

            //for (int i = 0; i< f_global.size(); i++) if (rank==0) std::cout << i <<" f global" << f_global[i]<< std::endl;
            //for (int i = 0; i< f.size(); i++) if (rank==0) std::cout << i <<" 0 f" << f[i]<< std::endl;
	    //MPI_Barrier(MPI_COMM_WORLD);
	    //for (int i = 0; i< f.size(); i++) if (rank==0) std::cout << i <<" 1 f local" << f[i]<< std::endl;
	    //MPI_Barrier(MPI_COMM_WORLD);
	    //std::exit(0);
            
            df.mv(u, newton_rhs);
            newton_rhs -= f;

            
            df_global.mv(u_global, newton_rhs_global);
            newton_rhs_global -= f_global;

            
            /*if(rank==0) std::cout << "*************************************************************** rechne Newton-Schritt " << i << std::endl;
            if(rank==0) std::cout << "globale Systemmatrix " << std::endl;
            for(int i=0; i< M_dune_global.N(); i++){
		for(int j=0; j< M_dune_global.M(); j++){
		if(rank==0) if(df_global.exists(i,j)){ std::cout << df_global[i][j]<< " ";}else{std::cout <<" ";}
		}if(rank==0)std::cout << std::endl;
	    }
	    
	    MPI_Barrier(MPI_COMM_WORLD);
            if(rank==0) std::cout << "lokale Systemmatrix " << std::endl;
	    for(int i=0; i< M_dune.N(); i++){
		for(int j=0; j< M_dune.M(); j++){
		if(rank==0) if(df.exists(i,j)){ std::cout << df[i][j]<< " ";}else{std::cout <<" ";}
		}if(rank==0)std::cout << std::endl;
	    }

	    MPI_Barrier(MPI_COMM_WORLD);

	    for(int i=0; i< M_dune.N(); i++){
		for(int j=0; j< M_dune.M(); j++){
		if(rank==1) if(df.exists(i,j)){ std::cout << df[i][j]<< " ";}else{std::cout <<" ";}
		}if(rank==1) std::cout << std::endl;
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
            
            if(rank==0) std::cout << "globale rechte seite " << std::endl;
            for (int i = 0; i< newton_rhs_global.size(); i++) if (rank==0) std::cout << i <<"rhs global" << newton_rhs_global[i]<< std::endl;
            if(rank==0) std::cout << "lokale rechte seite " << std::endl;
            for (int i = 0; i< newton_rhs.size(); i++) if (rank==0) std::cout << i <<"0 rhs local" << newton_rhs[i]<< std::endl;
	    MPI_Barrier(MPI_COMM_WORLD);
	    for (int i = 0; i< newton_rhs.size(); i++) if (rank==1) std::cout << i <<"1 rhs local" << newton_rhs[i]<< std::endl;
	    MPI_Barrier(MPI_COMM_WORLD);*/



///////////////////////////////////////////////////////////////////////////////////////

	//here we use parmg solver to solve the local system (df * u = newton_rhs) for the unknown u
  	auto &x=u;
  	x=0;
  	using MGSetup = Dune::ParMG::ParallelMultiGridSetup< FEBasis, MatrixType, VectorType >;
  	MGSetup mgSetup{*grid}; //0 grobgitter
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

//here we use dune-istl solver to solve the system (df_globl * u_global = newton_rhs_global)
	  
	  
	    //here we treat the boundary conditions
	    
	    //auto isDirichlet = [] (auto x) {return (x[0]<1e-8 or x[0]>0.9999 or x[1]<1e-8 or x[1]>0.9999);}; //two dimensional case

  	    auto isDirichlet = [] (auto x) {return (x[0] < -200.0 + 1e-8 or x[0] > 200.0-1e-8) ;};
  	    std::vector<char> dirichletNodes;
            interpolate(basis_global, dirichletNodes, isDirichlet);  

  	    for (size_t i=0; i<df_global.N(); i++){
    		if (dirichletNodes[i]){
      			auto cIt = df_global[i].begin();
      			auto cEndIt = df_global[i].end();
      			for(; cIt!=cEndIt; ++cIt){
        			*cIt = (i==cIt.index()) ? 1.0 : 0.0; 
      			}
    		}
  	    }
	    
	    //df_global[0][0]=1;
            //df_global[0][1]=0;

            //df_global[df_global.N()-1][df_global.M()-1]=1;
            //df_global[df_global.N()-1][df_global.M()-2]=0;
            
            newton_rhs_global[0] =0;
            newton_rhs_global[u_global.size()-1] =1;
          
          
            Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df_global);
            Dune::SeqILU<MatrixType,VectorType,VectorType> preconditioner(df_global, 0, 1.0);	  
            Dune::CGSolver<VectorType> cg_global(linearOperator,
                              preconditioner,
                              1e-16, 
                              5000,    
                              0);    
          
            Dune::InverseOperatorResult statistics_global ;
            newton_rhs_global[0]=0; newton_rhs_global[rhs_global.size()-1]=1;
            
            cg_global.apply(u_global, newton_rhs_global , statistics_global ); //rhs is not const



	    //for (int k=0; k<u_global.size(); k++)  if (rank==0) std::cout << "Newton " << i << " Index " << k<< " global u " << u_global[k] << " " << std::endl;	
	    //for (int k=0; k<u.size(); k++)  if (rank==0) std::cout << "rank 0 Newton " << i << " Index " << k<< " local u " << u[k] << " " << std::endl;
	    //MPI_Barrier(MPI_COMM_WORLD);
	    //for (int k=0; k<u.size(); k++)  if (rank==num_pro-1) std::cout << "rank 2 Newton " << i << " Index " << k<< " local u " << u[k] << " " << std::endl;
	    //MPI_Barrier(MPI_COMM_WORLD);

            evaluate_f(f_global, u_global, rhs_global, M_dune_global, A_dune_global, *w_global);   
            evaluate_f(f, u, rhs, M_dune, A_dune, *w); 
                     
                     
            using Norm =  EnergyNorm<MatrixType,VectorType>;
	    auto norm_global = Norm(A_dune_global); 
	    auto norm = Norm(A_dune);   
	    
	    
	    auto parallel_energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(A_dune, restrictToMaster, gridView.grid().comm());
	    
            if (rank==0) std::cout << i << " residuumsnorm von f_global(u) infinity " << norm_global(f_global) << std::endl;  	                  
            if (rank==0) std::cout << i << " residuumsnorm von f_global(u) energy " << f_global.infinity_norm() << std::endl;  
            std::cout << i << " residuumsnorm von f(u) " << parallel_energyNorm(f) << std::endl;  
            
            //if(rank==2) for(int i=0; i<rhs.size(); i++){std::cout << "0 *****  " << u[i] << std::endl;}	
            //if(f_global.infinity_norm()<1e-10){   std::cout << "genauigkeit erreicht global " << i << std::endl;      break;} 	  
            //if (rank==1) if(f.infinity_norm()<1e-10){   std::cout << "genauigkeit erreicht local" << i << std::endl;      break;} 	  

	    

	}
	
	
	


	/*for(int i=0; i< M_dune.N(); i++){
		for(int j=0; j< M_dune.M(); j++){
		if(M_dune.exists(i,j)){ std::cout << M_dune[i][j]<< " ";}else{std::cout <<" ";}
		}std::cout << std::endl;
	}*/

	//for(int i=0; i<u.size(); i++){std::cout << "***** Loesung " << u[i] << " exact "<< exact[i] << std::endl;}

	/*std::shared_ptr<TransferOperatorAssembler<GridType>> transfer;
	std::shared_ptr<TransferOperatorAssembler<GridType>> transfer_global;
	std::shared_ptr<std::vector<MatrixType*>> transferMatrix;
	std::shared_ptr<std::vector<MatrixType*>> transferMatrix_global;  
	
	transfer = std::make_shared<TransferOperatorAssembler<GridType>>(*grid);
	transfer_global = std::make_shared<TransferOperatorAssembler<GridType>>(*grid_global);
	transferMatrix = std::make_shared<std::vector<MatrixType*>>();
	transferMatrix_global = std::make_shared<std::vector<MatrixType*>>();
	for (int i=0; i< n_levels-1; i++){
	      transferMatrix->push_back(new MatrixType()); 
	      transferMatrix_global->push_back(new MatrixType()); 
	}
	transfer->assembleMatrixHierarchy<MatrixType>(*transferMatrix);
	transfer_global->assembleMatrixHierarchy<MatrixType>(*transferMatrix_global);
	MPI_Barrier(MPI_COMM_WORLD);
	
	std::shared_ptr<std::vector<MatrixType*>> vecvec = transferMatrix;
	std::shared_ptr<std::vector<MatrixType*>> vecvec_global = transferMatrix_global;

	
	MatrixType interpolate_matrix= *(vecvec->at(0));
	MatrixType restrict_matrix   = *(vecvec->at(0));

	MatrixType interpolate_matrix_global= *(vecvec_global->at(0));
	MatrixType restrict_matrix_global   = *(vecvec_global->at(0));
	
	//std::cout << "nach der ausgabe 7" << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	
	//std::cout << interpolate_matrix.N() <<" "<< interpolate_matrix.M() <<" "<< restrict_matrix.M() <<" "<< restrict_matrix.N() << std::endl;
	//std::cout << interpolate_matrix_global.N() <<" "<< interpolate_matrix_global.M() <<" "<< restrict_matrix_global.M() <<" "<< restrict_matrix_global.N() << std::endl;
	//std::cout << fe_basis[n_levels-2]->size() << " " << fe_basis[n_levels-1]->size() << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	//std::exit(0);
	for (int i=0; i< restrict_matrix.N(); i++){
		for (int j=0; j< restrict_matrix.M(); j++){
			if(restrict_matrix.exists(i,j)){	
		  		if (restrict_matrix[i][j]==0.5 ) restrict_matrix[i][j]=0;
			}

	      	}
	}
	
	for (int i=0; i< restrict_matrix_global.N(); i++){
		for (int j=0; j< restrict_matrix_global.M(); j++){
			if(restrict_matrix_global.exists(i,j)){	
		  		if (restrict_matrix_global[i][j]==0.5 ) restrict_matrix_global[i][j]=0;
			}

	      	}
	}
	


	restrict_matrix.mtv(u, coarse);
	restrict_matrix_global.mtv(u_global, coarse_global);
	interpolate_matrix.mtv(u, coarse);
	interpolate_matrix_global.mtv(u_global, coarse_global);
	
	interpolate_matrix.mv(coarse, u);
	interpolate_matrix_global.mv(coarse_global, u_global);*/

	//for (int i=0; i<coarse_global.size(); i++)  if (rank==0) std::cout <<i << " global u " << coarse_global[i] << " " << std::endl;
	//for (int i=0; i<coarse.size(); i++)  if (rank==0) std::cout <<i << " local u " << coarse[i] << " " << std::endl;		
	//for (int i=0; i<u_global.size(); i++)  if (rank==0) std::cout <<i << " interpolated global u " << u_global[i] << " " << std::endl;	
	//for (int i=0; i<u.size(); i++)  if (rank==0) std::cout <<i << " interpolated local u " << u[i] << " " << std::endl;
	
	MPI_Finalize();
}




