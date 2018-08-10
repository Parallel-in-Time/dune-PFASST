
#include ”config.h”
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

int main(int argc, char∗ argv[])
{
    MPIHelper::instance(argc,argv);
    auto world_comm = Dune::MPIHelper::getCollectiveCommunication();

    typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;
    typedef BlockVector<FieldVector<double,1> > VectorType;
    MatrixType A;
    VectorType b;

    // read matrix on rank 0
    if (world_comm.rank()==0 or true)
    {
        loadMatrixMarket(A, ”poisson−matrix.mm”);
        loadMatrixMarket(b, ”poisson−rhs.mm”);
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

    Communication∗ comm_redist;
    MatrixType parallel_A;
    typedef Dune::Amg::MatrixGraph<MatrixType> MatrixGraph;
    Dune::RedistributeInformation<Communication> rinfo;


    bool hasDofs = Dune::graphRepartition(MatrixGraph(A), comm,
                static_cast<int>(world_comm.size()),
                comm_redist,
                rinfo.getInterface (),
                true); // verbose

    rinfo.setSetup();
    redistributeMatrix(A, parallel_A , comm, ∗comm_redist, rinfo);

    VectorType parallel_b(parallel_A .N());
    VectorType parallel_x(parallel_A .M());
    rinfo.redistribute(b, parallel_b );

    // { problem setup end }


    // { solve begin }
    //Now we have the complete linear system and can solve it in parallel . To do this we need to set up the
    // parallel solver components and call the apply method of the solver . Below we use the conjugate gradient
    // method preconditioned with hybrid SSOR:

    if (hasDofs) // if hasDofs is false we do not compute.
    {
    // the index set has changed. Rebuild the remote information
    comm_redist−>remoteIndices().rebuild<false>();
    typedef Dune::SeqSSOR<MatrixType,VectorType,VectorType> Prec;
    typedef Dune::BlockPreconditioner<VectorType,VectorType, Communication,Prec> ParPrec; // type of parallel preconditioner
    typedef Dune::OverlappingSchwarzScalarProduct<VectorType,Communication> ScalarProduct; // type of parallel scalar product
    typedef Dune::OverlappingSchwarzOperator<MatrixType,VectorType, VectorType,Communication> Operator; // type of parallel linear operator

    ScalarProduct sp(∗comm_redist);

    Operator op(parallel_A, ∗comm redist);
    Prec prec(parallel_A , 1, 1.0);
    ParPrec pprec(prec, ∗comm_redist);

    // Object storing some statistics about the solving process
    Dune::InverseOperatorResult statistics ;
    Dune::CGSolver<VectorType> cg(op, // linear operator
                            sp,// scalar product
                            pprec,// parallel preconditioner
                            10e−8,// desired residual reduction factor
                            80,// maximum number of iterations
                            world_comm.rank()==0?2:0);// verbosity of the solver


    VectorType parallel_x(parallel_A.M());
    cg.apply( parallel_x , parallel_b , statistics );

    }

    // If you really need to then you can also gather all the information on the master process afterwards:

    VectorType x(A.M());
    rinfo .redistributeBackward(x, parallel_x );
    // Please note that you need to compile your program with DUNE’s ParMETIS flags.
    // { solve_end }
    return 0;
}
