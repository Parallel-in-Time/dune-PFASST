
/*#include <config.h>
#include <iostream>
#include <vector>
//


#include <dune/grid/yaspgrid.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/common/densematrix.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>

//#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
//#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
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
#include<iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/matrixmarket.hh>
#include<dune/istl/paamg/pinfo.hh>
#include<dune/istl/matrixredistribute.hh>
#include<dune/istl/paamg/graph.hh>

using namespace Dune;*/


// { problem setup begin }


/*using std::shared_ptr;
typedef Dune::YaspGrid<1> GridType; 
typedef GridType::LeafGridView GridView;
//using BaseFunction = Functions::PQkNodalBasis<GridView, 1>;
using BaseFunction = Functions::PQkNodalBasis<GridView, 1>;
const size_t nelements = 100;*/





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
typedef Dune::YaspGrid<1> GridType; 
typedef GridType::LeafGridView GridView;
using BaseFunction = Dune::Functions::PQkNodalBasis<GridView, 1>;
const size_t nelements = 100;






#define FIELD_DIM 1 // defines if the FieldVector is a scalar or a vector
#define SPACE_DIM 1 // defines if the problem is in 1D, 2D or 3D space







typedef double num_type; // to switch between single and double precision floating point numbers
using DuneVectorType = Dune::BlockVector<Dune::FieldVector<num_type, FIELD_DIM>>;
using DuneMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<num_type, FIELD_DIM, FIELD_DIM>>;


using namespace Dune;

int main(int argc, char* argv[])
{
    MPIHelper::instance(argc,argv);
    auto world_comm = Dune::MPIHelper::getCollectiveCommunication();

    typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;
    typedef BlockVector<FieldVector<double,1> > VectorType;
    MatrixType A;
    VectorType b;
    Dune::FieldVector<num_type,SPACE_DIM> x_start,x_end;
    std::array<int,SPACE_DIM> num_elements;
	

	for(int i = 0; i<SPACE_DIM; ++i)
	{
		x_start[i] = 0.0;
		x_end[i] = 1.0;
		num_elements[i] = 10;
	}

    // read matrix on rank 0
    //if (world_comm.rank()==0)
    
 		//std::cout<<"process with rank : "<<dune_comm.rank()<<" is here "<<std::endl;
		

	Dune::FieldVector<double,1> h = {1}; 
	std::array<int,1> n; 
	std::fill(n.begin(), n.end(), 10);
	std::shared_ptr<GridType> grid;
#if HAVE_MPI
    	grid = std::make_shared<GridType>(h, n, std::bitset<1>{0ULL}, 1, MPI_COMM_SELF); //mit overlapping
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



		if(world_comm.rank()==0){
		const int size = basis->size();  		
		int nnz = 3*size - 2;//4 + 3*(size-2)
		A.setBuildMode(DuneMatrixType::random); 
		//	set the size of the matrix
		A.setSize (size, size, nnz); 		
		// initially set row size for each row
  		A.setrowsize(0,2);
  		A.setrowsize(size-1,2);
  		for(int i=1;i<size-1; ++i)
      		A.setrowsize(i,3);
  		// finalize row setup phase
  		A.endrowsizes();
  		// add column entries to rows
  		A.addindex(0,0);
		A.addindex(0,1);
  		A.addindex(size-1,size-1);
		A.addindex(size-1,size-2);
   		for(int i=1;i<size-1; ++i)
   		{
      		A.addindex(i,i-1);
      		A.addindex(i,i);
      		A.addindex(i,i+1);
   		}
  		// finalize column setup phase
  		A.endindices();
  		// set entries using the random access operator
  		A[0][0] = num_elements[0];
		A[0][1] = -1*num_elements[0];
  		A[size-1][size-1] = num_elements[0];
		A[size-1][size-2] = -1*num_elements[0];
  		for(int i=1;i<size-1; ++i)
  		{
      		A[i][i-1] = -1.0*num_elements[0];
      		A[i][i]   = 2.0*num_elements[0];
      		A[i][i+1] = -1.0*num_elements[0];
  		}

		// filling the rhs with some values
		b.resize(num_elements[0]);
		for(int i=0;i <A.N(); ++i)
		{
			b[i] = i;
			std::cout << b[i]<<std::endl;
		}
		
		// displaying the matrix
	        		std::cout<<size << " " << std::endl;
		if(world_comm.rank()==0)
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(stiffness.exists(i,j))
          		{
	        		std::cout<<stiffness[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}
		}
    

    // Next we need the corresponding parallel index set . Note that this has to contain only mappings
    // for indices that might also be known to other processes. For matrix A this is not the case and
    // we can leave them empty for now:

    typedef std::size_t GlobalId; // The type for the global index
    typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;

    Communication comm(world_comm);
    // No need to add any indices to comm.indexSet()
    // Now we can create a parallel representation of the matrix and use PT−Scotch or ParMETIS
    // to distribute the sequential matrix to all processes . (See this post on how to install and use PT−Scotch.)

    Communication* comm_redist;
    MatrixType parallel_A;
    typedef Dune::Amg::MatrixGraph<MatrixType> MatrixGraph;
    Dune::RedistributeInformation<Communication> rinfo;

std::cout<<"Now Displaying the matrix part from process "<<std::endl;
    bool hasDofs = Dune::graphRepartition(MatrixGraph(A), comm,
                static_cast<int>(world_comm.size()),
                comm_redist,
                rinfo.getInterface(),
                true); // verbose
	std::cout<<"super "<<std::endl;
    rinfo.setSetup();		std::cout<<"Now Displaying the matrix part from process "<<std::endl;

    redistributeMatrix(A, parallel_A , comm, *comm_redist, rinfo);


for(int K=0; K<world_comm.size(); ++K){
	if(world_comm.rank() == K)
	{
		std::cout<<"Now Displaying the matrix part from process "<<world_comm.rank()<<std::endl;
		const int size = parallel_A.N();
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(parallel_A.exists(i,j))
          		{
	        		std::cout<<parallel_A[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}
	}	
	MPI_Barrier(world_comm);
	}//std::exit(0);



    VectorType parallel_b(parallel_A.N());
    VectorType parallel_x(parallel_A.M());
    rinfo.redistribute(b, parallel_b );

    // { problem setup end }


    // { solve begin }
    //Now we have the complete linear system and can solve it in parallel . To do this we need to set up the
    // parallel solver components and call the apply method of the solver . Below we use the conjugate gradient
    // method preconditioned with hybrid SSOR:

    if(hasDofs) // if hasDofs is false we do not compute.
    {
    // the index set has changed. Rebuild the remote information
    comm_redist->remoteIndices().rebuild<false>();

for(int K=0; K<world_comm.size(); ++K){
	if(world_comm.rank() == K)
	{
		std::cout<<"Now Displaying the matrix part from process "<<world_comm.rank()<<std::endl;
		const int size = parallel_A.N();
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(parallel_A.exists(i,j))
          		{
	        		std::cout<<parallel_A[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}
	}	
	MPI_Barrier(world_comm);
	}//std::exit(0);


    typedef Dune::SeqSSOR<MatrixType,VectorType,VectorType> Prec;
    typedef Dune::BlockPreconditioner<VectorType,VectorType, Communication,Prec> ParPrec; // type of parallel preconditioner
    typedef Dune::OverlappingSchwarzScalarProduct<VectorType,Communication> ScalarProduct; // type of parallel scalar product
    typedef Dune::OverlappingSchwarzOperator<MatrixType,VectorType, VectorType,Communication> Operator; // type of parallel linear operator

    ScalarProduct sp(*comm_redist);

    Operator op(parallel_A, *comm_redist);
    Prec prec(parallel_A , 1, 1.0);
    ParPrec pprec(prec, *comm_redist);

    // Object storing some statistics about the solving process
    Dune::InverseOperatorResult statistics ;
    Dune::CGSolver<VectorType> cg(op, // linear operator
                            sp,// scalar product
                            pprec,// parallel preconditioner
                            10e-8,// desired residual reduction factor
                            80,// maximum number of iterations
                            world_comm.rank()==0?2:0);// verbosity of the solver


    VectorType parallel_x(parallel_A.M());
    cg.apply( parallel_x , parallel_b , statistics );

    }

    // If you really need to then you can also gather all the information on the master process afterwards:

    /*VectorType x(A.M());
    rinfo.redistributeBackward(x, parallel_x );
		if(world_comm.rank()==0) for(int i=0;i <A.N(); ++i)
		{
			std::cout << x[i]<<std::endl;
		}*/
    // Please note that you need to compile your program with DUNE’s ParMETIS flags.
    // { solve_end }
    return 0;
}
