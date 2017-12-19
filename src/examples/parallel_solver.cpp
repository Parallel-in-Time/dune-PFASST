#include <config.h>
#include <iostream>
#include <vector>



#include <dune/grid/yaspgrid.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/common/densematrix.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>

//#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/typetree/utility.hh>

//#include <dune/fufem/assemblers/transferoperatorassembler.hh>



#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/geometry/quadraturerules.hh>

//#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>


#include <dune/istl/matrixmarket.hh>
#include <dune/istl/matrixredistribute.hh>

#include <dune/functions/functionspacebases/interpolate.hh>

#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh> 
#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>

//#include <dune/fufem/assemblers/istlbackend.hh>
//#include <dune/fufem/formatstring.hh>

#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/istlbackend.hh>
#include <dune/fufem/formatstring.hh>

#include <dune/istl/paamg/graph.hh> // Dune::Amg::MatrixGraph
#include <dune/istl/matrixredistribute.hh> // Dune::RedistributeInformation
#include <dune/istl/preconditioners.hh> // Dune::BlockPreconditioner, Dune::SeqSSOR
#include <dune/istl/solvers.hh> // Dune::CGSolver, Dune::RestartedGMResSolver
#include <dune/istl/schwarz.hh> // Dune::OverlappingSchwarzScalarProduct, Dune::OverlappingSchwarzOperator


#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication

//#include <dune/istl/ParallelRestrictedAdditiveSchwarz.hpp>

using std::shared_ptr;
typedef Dune::YaspGrid<2> GridType; 
typedef GridType::LeafGridView GridView;
using BaseFunction = Dune::Functions::PQkNodalBasis<GridView, 1>;
const size_t nelements = 8;


#ifndef PI
#define PI 3.1415926535897932385
#endif


int main(int argc, char** argv) {
    
  	MPI_Init(&argc, &argv);
        int my_rank, num_pro; MPI_Comm_rank(MPI_COMM_WORLD, &my_rank ); MPI_Comm_size(MPI_COMM_WORLD, &num_pro );


	Dune::FieldVector<double,2> h = {1,1}; std::array<int,2> n; std::fill(n.begin(), n.end(), nelements);
	std::shared_ptr<GridType> grid;
#if HAVE_MPI
    	grid = std::make_shared<GridType>(h, n, std::bitset<2>{0ULL}, 1, MPI_COMM_WORLD); //mit overlapping
#else
    	grid = std::make_shared<GridType>(h, n);
#endif    
        std::shared_ptr<BaseFunction> basis= std::make_shared<BaseFunction>(grid->leafGridView()); 


        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > mass, stiffness;
	auto sBackend = Dune::Fufem::istlMatrixBackend(stiffness);
        auto mBackend = Dune::Fufem::istlMatrixBackend(mass);
        using Assembler = Dune::Fufem::DuneFunctionsOperatorAssembler<BaseFunction, BaseFunction>;
        auto assembler = Assembler{*basis, *basis};
        using FiniteElement = std::decay_t<decltype(basis->localView().tree().finiteElement())>;

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
                

//rhs

        /*Dune::BlockVector<Dune::FieldVector<double,1> > initial(basis->size());
	const double t = 0;
        auto exact_solution = [t](const Dune::FieldVector<double,1>  &x){
            double solution=1.0;
            for(int i=0; i<1; i++){solution *= std::sin(PI * x[i]);}
            return solution * std::exp(-t * PI*PI);
        };

        interpolate(*basis, initial, exact_solution);*/

//rhs



	/*std::cout << "das ist jetzt das rhs " <<  std::endl;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );
	MPI_Barrier(MPI_COMM_WORLD);  	
	if(rank==0) for (size_t i = 0; i < rhs->get_data().size(); i++) {
          std::cout << "result " << rhs->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
	if(rank==1) for (size_t i = 0; i < rhs->get_data().size(); i++) {
          std::cout << "result " << rhs->data()[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD); 
	//std::exit(0);*/

        Dune::BlockVector<Dune::FieldVector<double,1> > u(basis->size());
        Dune::BlockVector<Dune::FieldVector<double,1> > M_rhs_dune ;
        //M_rhs_dune.resize(initial.size());
	//mass.mv(initial, M_rhs_dune);



       
	double dt=0.1;
        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dtA_dune = 	Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(stiffness);
        M_dtA_dune *= dt;
        M_dtA_dune += mass;





        std::cout << "das ist jetzt die matrix " <<  std::endl;	
	if(my_rank==0) 
    	Dune::printmatrix(std::cout, M_dtA_dune, "A", "");
	MPI_Barrier(MPI_COMM_WORLD);

	/*for (size_t i = 0; i < M_dtA_dune.M(); i++) {
			for (size_t j = 0; j < M_dtA_dune.N(); j++){
          			if (M_dtA_dune.exists(i,j)) {std::cout <<  M_dtA_dune[i][j] << " ";}else{std::cout  << 0 << " ";} }std::cout << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);*/
	/*if(rank==1)for (size_t i = 0; i < M_dtA_dune.M(); i++) {
			for (size_t j = 0; j < M_dtA_dune.N(); j++){
          			if (M_dtA_dune.exists(i,j)) {std::cout <<  M_dtA_dune[i][j] << " ";}else{std::cout  << 0 << " ";} }std::cout << std::endl;
        } 
	MPI_Barrier(MPI_COMM_WORLD); 
	if(rank==2)for (size_t i = 0; i < M_dtA_dune.M(); i++) {
			for (size_t j = 0; j < M_dtA_dune.N(); j++){
          			if (M_dtA_dune.exists(i,j)) {std::cout <<  M_dtA_dune[i][j] << " ";}else{std::cout  << 0 << " ";} }std::cout << std::endl;
        } 
	MPI_Barrier(MPI_COMM_WORLD); */
	std::exit(0);

	//MPI_Comm comm_x=MPI_COMM_WORLD; 
 	
	//MPI_Comm dune_comm = comm_x;//MPI_COMM_WORLD;
		
	using DuneVectorType = Dune::BlockVector<Dune::FieldVector<double, 1>>;
	using DuneMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;

	//using DuneCommunication = Dune::OwnerOverlapCopyCommunication<Dune::MatrixIndexSet, Dune::MatrixIndexSet>; //std::size_t>; decltype





	//auto localIndexSet = basis->localIndexSet();
	//auto localIndexSet2 = basis->localIndexSet();


	using IdSet = typename GridType::LocalIdSet;
	auto& idSet = grid->leafGridView().grid().localIdSet();

	using gIdSet = typename GridType::GlobalIdSet;
	auto& gidSet = grid->leafGridView().grid().globalIdSet();


	//GridType::LocalIdSet l = basis->GridType::IdSet::localIdSet;
	//auto g = basis->globalIdSet;

	using DuneCommunication = Dune::OwnerOverlapCopyCommunication<std::size_t>; //decltype(gidSet), decltype(idSet)>; //


	DuneCommunication DuneComm(MPI_COMM_WORLD);//dune_comm);
	DuneCommunication *comm_redist;

	Dune::MatrixIndexSet occupationPattern;

	occupationPattern.import(M_dtA_dune);

	//Dune::getOccupationPattern(basis, occupationPattern);

	//Dune::IndexInfoFromGrid<Dune::MatrixIndexSet, Dune::MatrixIndexSet> ifg; //(occupationPattern, localIndexSet);
	//Dune::IndexInfoFromGrid<decltype(gidSet), decltype(idSet)> ifg = Dune::IndexInfoFromGrid<decltype(gidSet), decltype(idSet)>(); // idSet, gidSet
	//ifg.Dune::IndexInfoFromGrid<decltype(gidSet), decltype(idSet)>::addLocalIndex(std::tuple<decltype(gidSet), decltype(idSet),int>(gidSet, idSet, 4));
	//ifg.Dune::IndexInfoFromGrid<decltype(gidSet), decltype(idSet)>::addRemoteIndex(std::tuple<int, decltype(gidSet), int>(my_rank, gidSet, 4));
	//ifg.addRemoteIndex(std::tuple<decltype(gidSet), decltype(idSet),int>(gidSet, idSet, 0);
	//DuneCommunication comm_new(ifg , MPI_COMM_WORLD);

	/*Dune::OwnerOverlapCopyCommunication< GlobalIdType, LocalIdType >::OwnerOverlapCopyCommunication	( 	const IndexInfoFromGrid< GlobalIdType, LocalIdType > &  	indexinfo,
		MPI_Comm  	comm_,
		SolverCategory::Category  	cat_ = SolverCategory::overlapping,
		bool  	freecomm_ = false 
	) */	


	using DuneMatrixGraph = Dune::Amg::MatrixGraph<DuneMatrixType>;
	Dune::RedistributeInformation<DuneCommunication> dune_rinfo;

	bool hasDofs = Dune::graphRepartition(DuneMatrixGraph(M_dtA_dune), DuneComm, static_cast<int>(num_pro),comm_redist, dune_rinfo.getInterface(), true);

	dune_rinfo.setSetup();
	//redistributeMatrix(M_dtA_dune, parallel_A, DuneComm, *comm_redist, dune_rinfo);



	//comm_redist->remoteIndices().rebuild<false>();
	using Seq_Preconditioner = Dune::SeqSSOR<DuneMatrixType, DuneVectorType, DuneVectorType>;			
	using Par_Preconditioner = Dune::BlockPreconditioner<DuneVectorType, DuneVectorType, DuneCommunication, Seq_Preconditioner>;
	using Par_ScalarProduct = Dune::OverlappingSchwarzScalarProduct<DuneVectorType, DuneCommunication>;
	using Par_LinearOperator = Dune::OverlappingSchwarzOperator<DuneMatrixType, DuneVectorType, DuneVectorType, DuneCommunication>;


	Par_ScalarProduct parallel_sp(*comm_redist);



	Par_LinearOperator parallel_linearoperator(M_dtA_dune, *comm_redist);

	Seq_Preconditioner seq_precon(M_dtA_dune, 1, 1.0); 
	Par_Preconditioner parallel_precon(seq_precon, *comm_redist);

	Dune::InverseOperatorResult statistics;



	Dune::RestartedGMResSolver<DuneVectorType> GMRES(parallel_linearoperator, parallel_sp,
														 parallel_precon,
														 1e-16,
														 200,
														 200,
														  my_rank ==0? 0:0);
		

	GMRES.apply(u, M_rhs_dune, statistics);

	MPI_Barrier(MPI_COMM_WORLD);  	
	if(my_rank==0) for (size_t i = 0; i < u.size(); i++) {
          std::cout << "nach lgs " << u[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==1) for (size_t i = 0; i < u.size(); i++) {
          std::cout << "nach lgs " << u[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank==2) for (size_t i = 0; i < u.size(); i++) {
          std::cout << "nach lgs " << u[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
	//std::exit(0);







   

	MPI_Finalize();
}


