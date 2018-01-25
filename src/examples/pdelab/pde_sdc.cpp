#include "config.h"
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>


// dune-geometry includes
//#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
//#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/newton/newton.hh>
//create mass + dt stiffness 
#include"nonlinearpoissonfem.hh"
//creates mass
#include"mass_operator.hh"
//constants
#ifndef PI
#define PI 3.1415926535897932385
#endif
const int DIM=1;
const int nelements=1000;




#include <dune/common/power.hh>
#include <dune/common/parametertree.hh>

#include <dune/istl/matrixmatrix.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istl/ovlpistlsolverbackend.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

//#include <dune/pdelab/finiteelementmap/pkqkfem.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>

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


#include "geometric_multigrid_components.hh"

      //using namespace Dune::PDELab::ISTL;

template<typename Number>
class Problem
{
public:
  typedef Number value_type;
  Problem () {}
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=DIM;
	const double t=0;	
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
        return solution * std::exp(-t * dim * PI*PI);
  }
};
template<typename Number>
class Problem2
{
public:
  typedef Number value_type;
  Problem2 () {}
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=DIM;
	const double t=0.2;	
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
        return solution * std::exp(-t * dim * PI*PI);
  }
};

int main(int argc, char** argv)
{

	Dune::MPIHelper&
      	helper = Dune::MPIHelper::instance(argc, argv);
    	//if(Dune::MPIHelper::isFake)
      		//std::cout<< "This is a sequential program." << std::endl;
    	//else
      		//std::cout << "Parallel code run on "    	<< helper.size() << " process(es)" << std::endl;

	double t1, t2; 
	t1 = MPI_Wtime();
	typedef double RF; 
	//setup the grid
        const int dim=DIM;
	const int degree =1;
        typedef Dune::YaspGrid<1> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<double,dim> h = {1};	      
	std::array<int,dim> n;
	std::fill(n.begin(), n.end(), nelements+1);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(h,n));

        gridp->refineOptions(false); // keep overlap in cells
        gridp->globalRefine(1);
        typedef Grid::LeafGridView GV;
	typedef Grid::LevelGridView GVl;
        GV gv=gridp->leafGridView();
        GVl gv_1=gridp->levelGridView(1);
        GVl gv_0=gridp->levelGridView(0);
	//typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
              //Dune::PDELab::QkLocalFiniteElementMap<GV, D, R, k >
	//typedef	Dune::PDELab::PkQkLocalFiniteElementMap<GV, DF, double> FEM;  
	typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF,double,  1> FEM;  
        FEM fem(gv);

  	//auto gt = Dune::GeometryTypes::quadrilateral;
  	//typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GVl::dimension> FEM;
  	//FEM fem(gt);


	std::cout << "hier" << std::endl;

	//typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM2;
        //FEM2 fem2(gv);


  	// Make grid function space
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	//typedef Dune::PDELab::GridFunctionSpace<GV,FEM2,CON,VBE> GFS2;
  	typedef Dune::PDELab::GridFunctionSpace<GVl,FEM,CON,VBE> GFSl;

  	GFS gfs(gv,fem);
  	//GFS2 gfs2(gv,fem2);


	std::cout << "hier 2" << std::endl;


	GFSl gfs_1(gv_1, fem);
	GFSl gfs_0(gv_0, fem);
  	using Zl = Dune::PDELab::Backend::Vector<GFSl,RF>;
	Zl v1(gfs_1);
	Zl v0(gfs_0);
	Zl coarse(gfs_1);


  	using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  	Z z(gfs); // mass times initial value
	Z initial(gfs); // initial value
	Z vgl(gfs); // analytic solution at the end point
	Z sol(gfs); // numeric solution
	Z neu(gfs);
	std::cout << "hier 3" << std::endl;

  	// Make a grid function out of it
  	typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  	ZDGF zdgf(gfs,z);

	std::cout << "hier 33" << std::endl;
  	Problem<RF> problem;
  	auto glambda = [&](const auto& x){return problem.g(x);};
  	auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);

	Problem2<RF> problem2;
  	auto glambda2 = [&](const auto& x){return problem2.g(x);};
  	auto g2 = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda2); 

	std::cout << "hier333" << std::endl;
  	// Fill the coefficient vector
  	Dune::PDELab::interpolate(g,gfs,initial);//z
  	Dune::PDELab::interpolate(g2,gfs,vgl); // vgl = analytic solution at the end point 

	//initial = z;

	//for(auto r =z.begin(); r !=z.end(); ++r){std::cout << "vector" << *r <<std::endl;}
 
		std::cout << "hier 4" << std::endl;

  	// make vector consistent NEW IN PARALLEL
  	Dune::PDELab::istl::ParallelHelper<GFS> grid_helper(gfs);
  	grid_helper.maskForeignDOFs(z);
  	Dune::PDELab::AddDataHandle<GFS,Z> adddh(gfs,z);
  	if (gfs.gridView().comm().size()>1){
    		gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
	}


  	// Assemble constraints
  	typedef typename GFS::template
    	ConstraintsContainer<RF>::Type CC;
  	CC cc;
  	
	// Make a local operator
  	typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;
  	LOP lop(problem);

  	// Make a global operator
  	typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  	MBE mbe((int)pow(1+2*degree,dim));
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GO;
  	GO go(gfs,cc,gfs,cc,lop,mbe);

	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

	typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType2;
	/*typedef MBE MatrixType2;
 	std::shared_ptr<TransferOperatorAssembler<Dune::YaspGrid<1>>> transfer;
	std::shared_ptr<std::vector<MatrixType2*>> transferMatrix;
	transfer = std::make_shared<TransferOperatorAssembler<Dune::YaspGrid<1>>>(*gridp);
	transferMatrix = std::make_shared<std::vector<MatrixType2*>>();
	for (int i=0; i< 1; i++){
		transferMatrix->push_back(new MatrixType2()); 
	}
	transfer->assembleMatrixHierarchy<MatrixType2>(*transferMatrix);
	typedef Dune::PDELab::NonoverlappingOperator<GFS, MatrixType2, Z, Z> NOO;
	NOO novlpOperator(gfs, (*transferMatrix)[0][0]);

	novlpOperator.apply(initial, neu);
	if(my_rank==1)for(auto r =initial.begin(); r !=initial.end(); ++r){std::cout << " initial " << *r <<std::endl;} 
	if(my_rank==1)for(auto r =neu.begin(); r !=neu.end(); ++r){std::cout << " neu     " << *r <<std::endl;} 

	Dune::PDELab::DiscreteGridFunction<GFSl, Zl> dgf(gfs_1, v1);
	Dune::PDELab::interpolate(dgf,gfs_0,v0);*/

	//typedef MBE MatrixType2;
	Dune::PDELab::interpolate(g,gfs_0,v0);
	Dune::PDELab::interpolate(g,gfs_1,coarse);
 	std::shared_ptr<TransferOperatorAssembler<Dune::YaspGrid<1>>> transfer;
	std::shared_ptr<std::vector<MatrixType2*>> transferMatrix;
	transfer = std::make_shared<TransferOperatorAssembler<Dune::YaspGrid<1>>>(*gridp);
	transferMatrix = std::make_shared<std::vector<MatrixType2*>>();
	for (int i=0; i< 1; i++){
		transferMatrix->push_back(new MatrixType2()); 
	}
	transfer->assembleMatrixHierarchy<MatrixType2>(*transferMatrix);

	/*for(int i=0; i<((*transferMatrix)[0][0]).M(); i++){
		for(int j=0; j<((*transferMatrix)[0][0]).N(); j++){ 
			if ((*transferMatrix)[0][0].exists(i,j)) {
				std::cout << (*transferMatrix)[0][0][i][j][0][0]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl;
	}*/

	typedef Dune::PDELab::NonoverlappingOperator<GFSl, MatrixType2, Zl, Zl> NOO;
	NOO novlpOperator(gfs_0, (*transferMatrix)[0][0]);

	(*transferMatrix)[0][0].mv(Dune::PDELab::Backend::native(v0), Dune::PDELab::Backend::native(v1));

	//novlpOperator.apply(v0, v1);
	for(auto r =v0.begin(); r !=v0.end(); ++r){std::cout << " initial " << *r <<std::endl;} 
	for(auto r =v1.begin(); r !=v1.end(); ++r){std::cout << " neu     " << *r <<std::endl;} 
	for(auto r =coarse.begin(); r !=coarse.end(); ++r){std::cout << " coarse     " << *r <<std::endl;} 


	//Dune::PDELab::DiscreteGridFunction<GFSl, Zl> dgf(gfs_1, v1);
	//Dune::PDELab::interpolate(dgf,gfs_0,v0);


	//MBE kk((*transferMatrix)[0]);
  	/*typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MatrixType2,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GO2;
  	GO2 go2(gfs,cc,gfs,cc,lop,(*transferMatrix)[0]);*/

	//ProlongationOperator<GFSl> pgo(gfs_0, gfs_1);
	//ProlongationOperator<GFS> prolong(gfs, gfs);
	//std::cout << "prolong erstellt" << std::endl;



  	typedef massFEM<Problem<RF>,FEM> LOPm;//////////////////////////////////
  	LOPm lopm(problem);

  	// Make a global operator
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOPm,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GOm;
  	GOm gom(gfs,cc,gfs,cc,lopm,mbe);

	gom.jacobian_apply(initial, z);




	//auto a = gom.getmat();

  	// make coefficent Vectors
  	using X = Dune::PDELab::Backend::Vector<GFS,double>;
  	X x(gfs,0.0);

  	// represent operator as a matrix
	typedef typename GO::template MatrixContainer<RF>::Type M;
  	M m(go);
  	//std::cout << m.patternStatistics() << std::endl;
  	m = 0.0;
  	go.jacobian(x,m);


  	// make coefficent Vectors
  	using Xm = Dune::PDELab::Backend::Vector<GFS,double>;
  	Xm xm(gfs,0.0);

  	// represent operator as a matrix
	//typedef typename GO::template MatrixContainer<RF>::Type Mm;
  	M mm(gom);
  	//std::cout << m.patternStatistics() << std::endl;
  	mm = 0.0;
  	gom.jacobian(xm,mm);


	mm*=0.002;
	//m += m;


	/*using 	Jacobian = Dune::PDELab::Backend::Matrix< MBE, X,X,double >;////////////////
	Jacobian j;
	go.jacobian(x, j);*/

	Dune::PDELab::Backend::native(m)[0][0][0][0] = 1;
	
	Dune::PDELab::Backend::native(m)[0][1][0][0] = 0;


	Dune::PDELab::Backend::native(m)[Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1] = 1;
	Dune::PDELab::Backend::native(m)[Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-2][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1] = 0;

	/*Dune::PDELab::Backend::native(m)[0][0];
	for(int i=0; i<Dune::PDELab::Backend::native(m).M(); i++){
		for(int j=0; j<Dune::PDELab::Backend::native(m).N(); j++){ 
			if (Dune::PDELab::Backend::native(m).exists(i,j)) {
				std::cout << Dune::PDELab::Backend::native(m)[i][j][0][0]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl;
	}*/



  	// Select a linear solver backend NEW IN PARALLEL
  	typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS; ///////
  	int verbose=0;
  	if (gfs.gridView().comm().rank()==0) verbose=1;
  	LS ls(gfs,100,verbose);


	for (int i=0; i<100; i++){
		ls.apply(m, sol, z, 0);
		//for(auto r =sol.begin(); r !=sol.end(); ++r){std::cout << "vector " << *r <<std::endl;}
		gom.jacobian_apply(sol, z);
		Dune::PDELab::Backend::native(z)[0][0]=0;
		Dune::PDELab::Backend::native(z)[Dune::PDELab::Backend::native(m).M()-1][0]=0;
		//for(auto r =z.begin(); r !=z.end(); ++r){std::cout << "new rhs " << *r <<std::endl;}
		//std::exit(0);
	}


           /*// restrict defect to CG subspace
        CGY cgd(p.M());
        p.mtv(d,cgd);


        // prolongate correction
        p.mv(cgv,v);
        dgmatrix.mmv(v,d);*/




	// M um+1 = M um + dt A um+1
	// (M + dt A) um+1 = M um 


	vgl -= sol; 
	for(auto r =vgl.begin(); r !=vgl.end(); ++r){std::cout << "vector " << *r <<std::endl;}
	//gom.applyscaleadd(-1, sol, vgl);
	std::cout << "norm" << Dune::PDELab::Backend::native(vgl).infinity_norm() << std::endl;

	for(int i=0; i<Dune::PDELab::Backend::native(m).N(); i++)
				std::cout << "initial " << Dune::PDELab::Backend::native(initial)[i][0] << " num " << Dune::PDELab::Backend::native(sol)[i][0] << " alg " << Dune::PDELab::Backend::native(vgl)[i][0] << std::endl;

	t2 = MPI_Wtime(); ///////////////
	printf( "Elapsed time is %f\n", t2 - t1 ); 


  	// Write VTK output file
  	Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
  	typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  	vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
                                         VTKF(zdgf,"fesol")));
  	vtkwriter.write("output",
        Dune::VTK::appendedraw);








	return 0;
}
