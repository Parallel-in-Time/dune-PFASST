/*#include <iostream>
#include <cmath>


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif


#include <dune/istl/bvector.hh> // BlockVector
#include <dune/istl/bcrsmatrix.hh> // BCRSMatrix

#include <dune/grid/yaspgrid.hh> // YaspGrid
#include <dune/functions/functionspacebases/pq1nodalbasis.hh> // PQ1NodalBasis
#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh> // DuneFunctionsOperatorAssembler
#include <dune/fufem/assemblers/istlbackend.hh> // istlMatrixBackend
#include <dune/fufem/assemblers/localassemblers/massassembler.hh> //MassAssembler
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh> //LaplaceAssembler
#include <dune/istl/paamg/graph.hh> // Dune::Amg::MatrixGraph
#include<dune/istl/matrixredistribute.hh> // Dune::RedistributeInformation


#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication

#include "mpi.h"




#define FIELD_DIM 1 // defines if the FieldVector is a scalar or a vector
#define SPACE_DIM 1 // defines if the problem is in 1D, 2D or 3D space



typedef double num_type; // to switch between single and double precision floating point numbers
using DuneVectorType = Dune::BlockVector<Dune::FieldVector<num_type, FIELD_DIM>>;
using DuneMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<num_type, FIELD_DIM, FIELD_DIM>>;


int main(int argc, char* argv[])
{

	Dune::MPIHelper::instance(argc, argv);
	auto dune_comm = Dune::MPIHelper::getCollectiveCommunication();

//	int world_rank=-1, world_size=-1;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	
//	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


	
	DuneMatrixType A,M;
	DuneVectorType b;

	Dune::FieldVector<num_type,SPACE_DIM> x_start,x_end;
	std::array<int,SPACE_DIM> num_elements;
	

	for(int i = 0; i<SPACE_DIM; ++i)
	{
		x_start[i] = 0.0;
		x_end[i] = 1.0;
		num_elements[i] = 32;
	}

	std::cout<<"DUNE rank is "<<dune_comm.rank()<<"\tof\t"<<dune_comm.size()<<std::endl;

//	std::cout<<"WORLD rank is "<<world_rank<<"\tof\t"<<world_size<<std::endl;

	//generate the FE-Matrix for the master process
	if(dune_comm.rank() == 0)
	{
		std::cout<<"process with rank : "<<dune_comm.rank()<<" is here "<<std::endl;
		


		const int size = num_elements[0] + 1;  		
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
		for(int i=0;i <M.N(); ++i)
		{
			b[i] = 1.0;
		}
		
		// displaying the matrix
		
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(A.exists(i,j))
          		{
	        		std::cout<<A[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}
	}

	using DuneCommunication = Dune::OwnerOverlapCopyCommunication<std::size_t>;
	DuneCommunication DuneComm(dune_comm); 
		
	DuneCommunication *comm_redist;
	DuneMatrixType parallel_A;
	using DuneMatrixGraph = Dune::Amg::MatrixGraph<DuneMatrixType>;
	Dune::RedistributeInformation<DuneCommunication> dune_rinfo;
	
	bool hasDofs = Dune::graphRepartition(DuneMatrixGraph(A), 
										  DuneComm,
										  static_cast<int>(dune_comm.size()),
										  comm_redist,
										  dune_rinfo.getInterface(),
										  true);

	dune_rinfo.setSetup();
	redistributeMatrix(A, parallel_A, DuneComm, *comm_redist, dune_rinfo);


	for(int K=0; K<dune_comm.size(); ++K){
	if(dune_comm.rank() == K)
	{
		std::cout<<"Now Displaying the matrix part from process "<<dune_comm.rank()<<std::endl;
		const int size = parallel_A.M();
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
	MPI_Barrier(dune_comm);
	}

//	MPI_Finalize();


}*/


/*#include <config.h>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/paamg/pinfo.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/matrixmarket.hh>
#include<dune/istl/matrixredistribute.hh>
#include<dune/istl/paamg/graph.hh>
#include<dune/common/parallel/mpihelper.hh>


int main(int argc, char** argv)
{

Dune::MPIHelper::instance(argc, argv);
auto world_comm = Dune::MPIHelper::getCollectiveCommunication();
const int BS=1; // block size, sparse scalar matrix
typedef Dune::FieldMatrix<double,BS,BS> MatrixBlock; // matrix block type
typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat; // sparse matrix type
typedef Dune::FieldVector<double,BS> VectorBlock; // vector block type      
typedef Dune::BlockVector<VectorBlock> BVector; // vector type
BCRSMat A;

// read matrix on rank 0
if(world_comm.rank()==0)
{}



typedef std::size_t GlobalId; // The type for the global index
typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;

Communication comm(world_comm);
Communication comm_world(Dune::MPIHelper::getCollectiveCommunication());
// No need to add any indices to comm.indexSet()



Communication* comm_redist;
BCRSMat parallel_A;
typedef Dune::Amg::MatrixGraph<BCRSMat> MatrixGraph;
Dune::RedistributeInformation<Communication> rinfo;
bool hasDofs = Dune::graphRepartition(MatrixGraph(A), comm,
                                      static_cast<int>(world_comm.size()),
                                      comm_redist,
                                      rinfo.getInterface(),
                                      true);
rinfo.setSetup();
redistributeMatrix(A, parallel_A, comm, comm_world, rinfo); //comm_world

BVector b(A.N());
BVector parallel_b(parallel_A.N());
BVector parallel_x(parallel_A.M());
b = 100.0;
rinfo.redistribute(b, parallel_b);



if(hasDofs) // if hasDofs is false we do not compute.
{
  // the index set has changed. Rebuild the remote information
  comm_redist->remoteIndices().rebuild<false>();
  typedef Dune::SeqSSOR<BCRSMat,BVector,BVector> Prec;
  //typedef Dune::BlockPreconditioner<VectorType,VectorType, Communication,Prec> ParPrec; // type of parallel preconditioner
  typedef Dune::BlockPreconditioner<BVector,BVector,Communication,Prec> ParPrec; // type of parallel preconditioner 
  typedef Dune::OverlappingSchwarzScalarProduct<BVector,Communication> ScalarProduct; // type of parallel scalar product
  typedef Dune::OverlappingSchwarzOperator<BCRSMat,BVector, BVector,Communication> Operator; // type of parallel linear operator

  ScalarProduct sp(*comm_redist);
  Operator op(parallel_A, *comm_redist);
  Prec prec(parallel_A, 1, 1.0);
  ParPrec pprec(prec, *comm_redist);
  Dune::InverseOperatorResult r;
  Dune::CGSolver<BVector> cg(op, sp, pprec, 10e-8, 80,
                             world_comm.rank()==0?2:0);
  BVector parallel_x(parallel_A.M());
  cg.apply(parallel_x, parallel_b, r);
}



BVector x(A.M());
ri.redistributeBackwards(x, parallel_x);



add_dune_parmetis_flags(<target-name>)




}*/
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*#include "config.h"

#define DEBUG_REPART

#include <dune/istl/matrixredistribute.hh>
#include <iostream>
#include <dune/istl/paamg/test/anisotropic.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixutils.hh>
#include <dune/istl/paamg/graph.hh>
#include <dune/istl/io.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/bigunsignedint.hh>

class MPIError {
public:

  MPIError(std::string s, int e) : errorstring(s), errorcode(e){}

  std::string errorstring;

  int errorcode;
};

void MPI_err_handler(MPI_Comm *, int *err_code, ...){
  char *err_string=new char[MPI_MAX_ERROR_STRING];
  int err_length;
  MPI_Error_string(*err_code, err_string, &err_length);
  std::string s(err_string, err_length);
  std::cerr << "An MPI Error occurred:"<<std::endl<<s<<std::endl;
  delete[] err_string;
  throw MPIError(s, *err_code);
}

template<int BS>
int testRepart(int N, int coarsenTarget)
{

  std::cout<<"==================================================="<<std::endl;
  std::cout<<"BS="<<BS<<" N="<<N<<" coarsenTarget="<<coarsenTarget<<std::endl;

  int procs, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  typedef Dune::FieldMatrix<double,BS,BS> MatrixBlock;
  typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
  typedef Dune::bigunsignedint<56> GlobalId;
  typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;
  int n;

  N/=BS;

  Communication comm(MPI_COMM_WORLD);

  BCRSMat mat = setupAnisotropic2d<BCRSMat>(N, comm.indexSet(), comm.communicator(), &n, 1);
  typedef typename Dune::Amg::MatrixGraph<BCRSMat> MatrixGraph;

  MatrixGraph graph(mat);
  Communication* coarseComm;

  comm.remoteIndices().template rebuild<false>();

  std::cout<<comm.communicator().rank()<<comm.indexSet()<<std::endl;

  Dune::RedistributeInformation<Communication> ri;
  Dune::graphRepartition(graph, comm, coarsenTarget,
                         coarseComm, ri.getInterface());

  std::cout<<coarseComm->communicator().rank()<<coarseComm->indexSet()<<std::endl;
  BCRSMat newMat;

  if(comm.communicator().rank()==0)
    std::cout<<"Original matrix"<<std::endl;
  comm.communicator().barrier();
  //printGlobalSparseMatrix(mat, comm, std::cout);
  const int size = mat.M();
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(mat.exists(i,j))
          		{
	        		std::cout<<mat[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}

  redistributeMatrix(mat, newMat, comm, *coarseComm, ri);
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout<<" groesse " << mat.M() << " " << newMat.M() <<std::endl;
  std::cout<<" new mat" <<std::endl;
  if(coarseComm->communicator().rank()==0){
  const int size = newMat.M();
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(newMat.exists(i,j))
          		{
	        		std::cout<<newMat[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}
  }
  std::cout<<" new mat" <<std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  std::cout<<comm.communicator().rank()<<": redist interface "<<ri.getInterface()<<std::endl;

  if(comm.communicator().rank()==0)
    std::cout<<"Redistributed matrix"<<std::endl;
  comm.communicator().barrier();
  if(coarseComm->communicator().size()>0)
    //printGlobalSparseMatrix(newMat, *coarseComm, std::cout);
  comm.communicator().barrier();

  // Check for symmetry
  int ret=0;
  typedef typename BCRSMat::ConstRowIterator RIter;
  for(RIter row=newMat.begin(), rend=newMat.end(); row != rend; ++row) {
    typedef typename BCRSMat::ConstColIterator CIter;
    for(CIter col=row->begin(), cend=row->end(); col!=cend; ++col)
    {
      if(col.index()<=row.index())
        try{
          newMat[col.index()][row.index()];
        }catch(Dune::ISTLError e) {
          std::cerr<<coarseComm->communicator().rank()<<": entry ("
                   <<col.index()<<","<<row.index()<<") missing!"<<std::endl;
          ret=1;

        }
      else
        break;
    }
  }

  //if(coarseComm->communicator().rank()==0)
  //Dune::printmatrix(std::cout, newMat, "redist", "row");
  delete coarseComm;
  return ret;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  MPI_Errhandler handler;
  MPI_Errhandler_create(MPI_err_handler, &handler);
  MPI_Errhandler_set(MPI_COMM_WORLD, handler);
  int procs;
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  int N=4*procs;

  int coarsenTarget=1;

  if(argc>1)
    N = atoi(argv[1]);
  if(argc>2)
    coarsenTarget = atoi(argv[2]);

  if(N<procs*2) {
    std::cerr<<"Problem size insufficient for process number"<<std::endl;
    return 1;
  }

  testRepart<1>(N,coarsenTarget);
  MPI_Finalize();
}*/

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*#include "config.h"
#ifdef TEST_AGGLO
#define UNKNOWNS 10
#endif
#include "anisotropic.hh"
#include <dune/common/timer.hh>
#include <dune/common/unused.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <string>

template<class T, class C>
class DoubleStepPreconditioner
  : public Dune::Preconditioner<typename T::domain_type, typename T::range_type>
{
public:
  typedef typename T::domain_type X;
  typedef typename T::range_type Y;

  enum {category = T::category};

  DoubleStepPreconditioner(T& preconditioner_, C& comm)
    : preconditioner(&preconditioner_), comm_(comm)
  {}

  virtual void pre (X& x, Y& b)
  {
    preconditioner->pre(x,b);
  }

  virtual void apply(X& v, const Y& d)
  {
    preconditioner->apply(v,d);
    comm_.copyOwnerToAll(v,v);
  }

  virtual void post (X& x)
  {
    preconditioner->post(x);
  }
private:
  T* preconditioner;
  C& comm_;
};


class MPIError {
public:

  MPIError(std::string s, int e) : errorstring(s), errorcode(e){}

  std::string errorstring;

  int errorcode;
};

void MPI_err_handler(MPI_Comm *comm, int *err_code, ...){
  DUNE_UNUSED_PARAMETER(comm);
  char *err_string=new char[MPI_MAX_ERROR_STRING];
  int err_length;
  MPI_Error_string(*err_code, err_string, &err_length);
  std::string s(err_string, err_length);
  std::cerr << "An MPI Error occurred:"<<std::endl<<s<<std::endl;
  delete[] err_string;
  throw MPIError(s, *err_code);
}

template<int BS>
void testAmg(int N, int coarsenTarget)
{
  std::cout<<"==================================================="<<std::endl;
  std::cout<<"BS="<<BS<<" N="<<N<<" coarsenTarget="<<coarsenTarget<<std::endl;

  int procs, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  typedef Dune::FieldMatrix<double,BS,BS> MatrixBlock;
  typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
  typedef Dune::FieldVector<double,BS> VectorBlock;
  typedef Dune::BlockVector<VectorBlock> Vector;
  typedef int GlobalId;
  typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;
  typedef Dune::OverlappingSchwarzOperator<BCRSMat,Vector,Vector,Communication> Operator;
  int n;

  N/=BS;

  Communication comm(MPI_COMM_WORLD);

  BCRSMat mat = setupAnisotropic2d<BCRSMat>(N, comm.indexSet(), comm.communicator(), &n, 1);

  const BCRSMat& cmat = mat;

  comm.remoteIndices().template rebuild<false>();

  Vector b(cmat.N()), x(cmat.M());

  b=0;
  x=100;

  setBoundary(x, b, N, comm.indexSet());

  Vector b1=b, x1=x;

  if(N<=6) {
    std::ostringstream name;
    name<<rank<<": row";

    Dune::printmatrix(std::cout, cmat, "A", name.str().c_str());
    Dune::printvector(std::cout, x, "x", name.str().c_str());
    //Dune::printvector(std::cout, b, "b", name.str().c_str());
    //Dune::printvector(std::cout, b1, "b1", "row");
    //Dune::printvector(std::cout, x1, "x1", "row");
  }

  Dune::Timer watch;

  watch.reset();
  Operator fop(cmat, comm);

  typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<BCRSMat,Dune::Amg::FirstDiagonal> >
  Criterion;
  typedef Dune::SeqSSOR<BCRSMat,Vector,Vector> Smoother;
  //typedef Dune::SeqJac<BCRSMat,Vector,Vector> Smoother;
  //typedef Dune::SeqILU0<BCRSMat,Vector,Vector> Smoother;
  //typedef Dune::SeqILUn<BCRSMat,Vector,Vector> Smoother;
  typedef Dune::BlockPreconditioner<Vector,Vector,Communication,Smoother> ParSmoother;
  typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;

  Dune::OverlappingSchwarzScalarProduct<Vector,Communication> sp(comm);

  Dune::InverseOperatorResult r, r1;

  double buildtime;

  SmootherArgs smootherArgs;

  smootherArgs.iterations = 1;


  Criterion criterion(15,coarsenTarget);
  criterion.setDefaultValuesIsotropic(2);


  typedef Dune::Amg::AMG<Operator,Vector,ParSmoother,Communication> AMG;

  AMG amg(fop, criterion, smootherArgs, comm);

  buildtime = watch.elapsed();

  if(rank==0)
    std::cout<<"Building hierarchy took "<<buildtime<<" seconds"<<std::endl;

  Dune::CGSolver<Vector> amgCG(fop, sp, amg, 10e-8, 300, (rank==0) ? 2 : 0);
  watch.reset();

  amgCG.apply(x,b,r);
  amg.recalculateHierarchy();


  MPI_Barrier(MPI_COMM_WORLD);
  double solvetime = watch.elapsed();

  b=0;
  x=100;

  setBoundary(x, b, N, comm.indexSet());

  Dune::CGSolver<Vector> amgCG1(fop, sp, amg, 10e-8, 300, (rank==0) ? 2 : 0);
  amgCG1.apply(x,b,r);

  if(!r.converged && rank==0)
    std::cerr<<" AMG Cg solver did not converge!"<<std::endl;

  if(rank==0) {
    std::cout<<"AMG solving took "<<solvetime<<" seconds"<<std::endl;

    std::cout<<"AMG building took "<<(buildtime/r.elapsed*r.iterations)<<" iterations"<<std::endl;
    std::cout<<"AMG building together with slving took "<<buildtime+solvetime<<std::endl;
  }

}

template<int BSStart, int BSEnd, int BSStep=1>
struct AMGTester
{
  static void test(int N, int coarsenTarget)
  {
    testAmg<BSStart>(N, coarsenTarget);
    const int next = (BSStart+BSStep>BSEnd) ? BSEnd : BSStart+BSStep;
    AMGTester<next,BSEnd,BSStep>::test(N, coarsenTarget);
  }
}
;

template<int BSStart,int BSStep>
struct AMGTester<BSStart,BSStart,BSStep>
{
  static void test(int N, int coarsenTarget)
  {
    testAmg<BSStart>(N, coarsenTarget);
  }
};


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  MPI_Errhandler handler;
  MPI_Errhandler_create(MPI_err_handler, &handler);
  MPI_Errhandler_set(MPI_COMM_WORLD, handler);

  int N=100;

  int coarsenTarget=200;

  if(argc>1)
    N = atoi(argv[1]);

  if(argc>2)
    coarsenTarget = atoi(argv[2]);

#ifdef TEST_AGGLO
  N=UNKNOWNS;
#endif
  AMGTester<1,1>::test(N, coarsenTarget);
  //AMGTester<1,5>::test(N, coarsenTarget);
  //  AMGTester<10,10>::test(N, coarsenTarget);

  MPI_Finalize();
}*/

/*#include <iostream>
#include <cmath>
#include <config.h>

#include <dune/common/densematrix.hh>
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/matrixredistribute.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/typetree/utility.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh> 

#include <dune/fufem/assemblers/transferoperatorassembler.hh>


//#ifdef HAVE_CONFIG_H
//#endif

#include "assemble.hpp"
 
#include <dune/istl/bvector.hh> // BlockVector
#include <dune/istl/bcrsmatrix.hh> // BCRSMatrix

#include <dune/grid/yaspgrid.hh> // YaspGrid
#include <dune/functions/functionspacebases/pq1nodalbasis.hh> // PQ1NodalBasis
#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh> // DuneFunctionsOperatorAssembler
#include <dune/fufem/assemblers/istlbackend.hh> // istlMatrixBackend
#include <dune/fufem/assemblers/localassemblers/massassembler.hh> //MassAssembler
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh> //LaplaceAssembler
#include <dune/istl/paamg/graph.hh> // Dune::Amg::MatrixGraph
#include<dune/istl/matrixredistribute.hh> // Dune::RedistributeInformation


#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication

#include <dune/istl/schwarz.hh>

#include "mpi.h"




#define FIELD_DIM 1 // defines if the FieldVector is a scalar or a vector
#define SPACE_DIM 1 // defines if the problem is in 1D, 2D or 3D space

const int BASE_ORDER=1;
const int nelements=10;

typedef double num_type; // to switch between single and double precision floating point numbers
using DuneVectorType = Dune::BlockVector<Dune::FieldVector<num_type, FIELD_DIM>>;
using DuneMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<num_type, FIELD_DIM, FIELD_DIM>>;


int main(int argc, char* argv[])
{

	Dune::MPIHelper::instance(argc, argv);
	auto dune_comm = Dune::MPIHelper::getCollectiveCommunication();

//	int world_rank=-1, world_size=-1;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	
//	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


	
	DuneMatrixType A;
        DuneMatrixType M;
	DuneVectorType b;

	Dune::FieldVector<num_type,SPACE_DIM> x_start,x_end;
	std::array<int,SPACE_DIM> num_elements;
	

	for(int i = 0; i<SPACE_DIM; ++i)
	{
		x_start[i] = 0.0;
		x_end[i] = 1.0;
		num_elements[i] = 32;
	}

	std::cout<<"DUNE rank is "<<dune_comm.rank()<<"\tof\t"<<dune_comm.size()<<std::endl;

//	std::cout<<"WORLD rank is "<<world_rank<<"\tof\t"<<world_size<<std::endl;

	//generate the FE-Matrix for the master process
           typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
    typedef GridType::LevelGridView GridView;
    using BasisFunction = Dune::Functions::PQkNodalBasis<GridView, BASE_ORDER>;
     
    //typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; 
    
        std::shared_ptr<GridType> grid;
    

    int n_levels=1;

    std::vector<std::shared_ptr<BasisFunction> > fe_basis(n_levels); ; 

    
    Dune::FieldVector<double,1> hR = {200};
    Dune::FieldVector<double,1> hL = {-200};
    std::array<int,1> n;
    std::fill(n.begin(), n.end(), nelements); 	    
//#if HAVE_MPI
    //grid = std::make_shared<GridType>(hL, hR, n, std::bitset<1>{0ULL}, 1, MPI_COMM_SELF);
//#else
    grid = std::make_shared<GridType>(hL, hR, n);
//#endif
    for (int i=0; i<n_levels; i++){	      
	      //grid->globalRefine((bool) i);
	      auto view = grid->levelGridView(i);
              fe_basis[n_levels-i-1] = std::make_shared<BasisFunction>(grid->levelGridView(i)); //grid->levelGridView(i));//gridView);
    } 
        
        
	if(dune_comm.rank() == 0)
	{
                    assembleProblem(fe_basis[0], A, M);

	}

	using DuneCommunication = Dune::OwnerOverlapCopyCommunication<std::size_t>;
	DuneCommunication DuneComm(dune_comm); 
		
	DuneCommunication *comm_redist;
	DuneMatrixType parallel_A;
	using DuneMatrixGraph = Dune::Amg::MatrixGraph<DuneMatrixType>;
	Dune::RedistributeInformation<DuneCommunication> dune_rinfo;
			std::cout<< "vor " << dune_comm.size() <<std::endl;
	
	bool hasDofs = Dune::graphRepartition(DuneMatrixGraph(A), 
										  DuneComm,
										  static_cast<int>(dune_comm.size()),
										  comm_redist,
										  dune_rinfo.getInterface(),
										  true);
			std::cout<< "nach" << std::endl;

	dune_rinfo.setSetup();
	redistributeMatrix(A, parallel_A, DuneComm, *comm_redist, dune_rinfo);



    DuneVectorType parallel_b(parallel_A.N());
    DuneVectorType parallel_x(parallel_A.M());
    dune_rinfo.redistribute(b, parallel_b );

    // { problem setup end }


    // { solve begin }
    //Now we have the complete linear system and can solve it in parallel . To do this we need to set up the
    // parallel solver components and call the apply method of the solver . Below we use the conjugate gradient
    // method preconditioned with hybrid SSOR:

    if (hasDofs) // if hasDofs is false we do not compute.
    {
    // the index set has changed. Rebuild the remote information
    comm_redist->remoteIndices().rebuild<false>();
    typedef Dune::SeqSSOR<DuneMatrixType,DuneVectorType,DuneVectorType> Prec;
    typedef Dune::BlockPreconditioner<DuneVectorType,DuneVectorType, DuneCommunication,Prec> ParPrec; // type of parallel preconditioner
    typedef Dune::OverlappingSchwarzScalarProduct<DuneVectorType,DuneCommunication> ScalarProduct; // type of parallel scalar product
    typedef Dune::OverlappingSchwarzOperator<DuneMatrixType,DuneVectorType, DuneVectorType, DuneCommunication> Operator; // type of parallel linear operator

    ScalarProduct sp(*comm_redist);

    Operator op(parallel_A, *comm_redist);
    Prec prec(parallel_A , 1, 1.0);
    ParPrec pprec(prec, *comm_redist);

    // Object storing some statistics about the solving process
    Dune::InverseOperatorResult statistics ;
    Dune::CGSolver<DuneVectorType> cg(op, // linear operator
                            sp,// scalar product
                            pprec,// parallel preconditioner
                            10e-8,// desired residual reduction factor
                            80,// maximum number of iterations
                            dune_comm.rank()==0?2:0);// verbosity of the solver


    DuneVectorType parallel_x(parallel_A.M());
    cg.apply( parallel_x , parallel_b , statistics );

    }

    // If you really need to then you can also gather all the information on the master process afterwards:

    DuneVectorType x(A.M());
    dune_rinfo.redistributeBackward(x, parallel_x );


	for(int K=0; K<dune_comm.size(); ++K){
	if(dune_comm.rank() == K)
	{
		std::cout<<"Now Displaying the matrix part from process "<<dune_comm.rank()<<std::endl;
		const int size = parallel_A.M();
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
	MPI_Barrier(dune_comm);
	}

//	MPI_Finalize();


}*/
#include <config.h>
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

using namespace Dune;


// { problem setup begin }


#define FIELD_DIM 1 // defines if the FieldVector is a scalar or a vector
#define SPACE_DIM 1 // defines if the problem is in 1D, 2D or 3D space



typedef double num_type; // to switch between single and double precision floating point numbers
using DuneVectorType = Dune::BlockVector<Dune::FieldVector<num_type, FIELD_DIM>>;
using DuneMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<num_type, FIELD_DIM, FIELD_DIM>>;

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
		num_elements[i] = 5;
	}

    // read matrix on rank 0
    if (world_comm.rank()==0 or true)
    {
 		//std::cout<<"process with rank : "<<dune_comm.rank()<<" is here "<<std::endl;
		


		const int size = num_elements[0] + 1;  		
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
		
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(A.exists(i,j))
          		{
	        		std::cout<<A[i][j] << " ";
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


    bool hasDofs = Dune::graphRepartition(MatrixGraph(A), comm,
                static_cast<int>(world_comm.size()),
                comm_redist,
                rinfo.getInterface (),
                true); // verbose

    rinfo.setSetup();
    redistributeMatrix(A, parallel_A , comm, *comm_redist, rinfo);

    VectorType parallel_b(parallel_A .N());
    VectorType parallel_x(parallel_A .M());
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

    VectorType x(A.M());
    rinfo.redistributeBackward(x, parallel_x );
		if(world_comm.rank()==0) for(int i=0;i <A.N(); ++i)
		{
			std::cout << x[i]<<std::endl;
		}
    // Please note that you need to compile your program with DUNE’s ParMETIS flags.
    // { solve_end }
    return 0;
}
