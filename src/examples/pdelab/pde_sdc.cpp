#include "config.h"
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
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

#include"nonlinearpoissonfem.hh"
#include"mass_operator.hh"

//#include"nonlinearheatfem.hh"
#ifndef PI
#define PI 3.1415926535897932385
#endif

const int DIM=1;
const int nelements=10;
/*template<typename Number>
class Problem
{

public:
  typedef Number value_type;

  //! Constructor without arg sets nonlinear term to zero
  Problem () {}

  //! Dirichlet extension
  template<typename X>
  Number g (const X& x) const
  {
	double solution=1.0;
        solution *= std::sin(PI * x.size());
        return solution;
  }

};*/

template<typename Number>
class Problem
{
public:
  typedef Number value_type;


  //! Constructor without arg sets nonlinear term to zero
  Problem () {}


  //! Dirichlet extension
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=DIM;
	const double t=0;	
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
	//std::cout << solution;
        return solution * std::exp(-t * dim * PI*PI);
  }


};


template<typename Number>
class Problem2
{
public:
  typedef Number value_type;


  //! Constructor without arg sets nonlinear term to zero
  Problem2 () {}


  //! Dirichlet extension
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=DIM;
	const double t=0.2;	
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
	//std::cout << solution;
        return solution * std::exp(-t * dim * PI*PI);
  }


};

int main(int argc, char** argv)
{

	Dune::MPIHelper&
      	helper = Dune::MPIHelper::instance(argc, argv);
    	if(Dune::MPIHelper::isFake)
      		std::cout<< "This is a sequential program." << std::endl;
    	else
      		std::cout << "Parallel code run on "
                	<< helper.size() << " process(es)" << std::endl;

	double t1, t2; 
	t1 = MPI_Wtime();


        const int dim=DIM;
	const int degree =1;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        /*Dune::FieldVector<DF,dim> L;
        L[0] = 1;
        L[1] = 1;
        std::array<int,dim> N;
        N[0] = 16;//16;
        N[1] = 16;
        std::bitset<dim> B(false);
        int overlap=1;
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication()));*/
	    Dune::FieldVector<double,dim> h = {1};
	    
	      
	    std::array<int,dim> n;
	    std::fill(n.begin(), n.end(), nelements);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(h,n));

        gridp->refineOptions(false); // keep overlap in cells
        //gridp->globalRefine(1);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        FEM fem(gv);
        //driver(gv,fem,ptree);

  	// Make grid function space
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON; // NEW IN PARALLEL
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	GFS gfs(gv,fem);


	typedef double RF;                   // type for computations
  	// A coefficient vector
  	using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  	Z z(gfs); // initial value
	Z initial(gfs);
	Z vgl(gfs);
	Z sol(gfs);

  	// Make a grid function out of it
  	typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  	ZDGF zdgf(gfs,z);

  	//Problem<RF> problem();
	//auto glambda = [&](const auto& x){return problem.g(x);};
  	//auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);
	//RF eta = 1.0;
  	Problem<RF> problem;
  	auto glambda = [&](const auto& x){return problem.g(x);};
  	auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);

	Problem2<RF> problem2;
  	auto glambda2 = [&](const auto& x){return problem2.g(x);};
  	auto g2 = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda2);

  	// Fill the coefficient vector
  	Dune::PDELab::interpolate(g,gfs,z);
  	Dune::PDELab::interpolate(g2,gfs,vgl);

	initial = z;
	//for(auto r =z.begin(); r !=z.end(); ++r){std::cout << "vector" << *r <<std::endl;}
 
	//z.data(0);
	//for(int i=0; i<z.size(); i++){}
	Dune::PDELab::Backend::native(z)[0];

  	// make vector consistent NEW IN PARALLEL
  	Dune::PDELab::istl::ParallelHelper<GFS> grid_helper(gfs);
  	grid_helper.maskForeignDOFs(z);
  	Dune::PDELab::AddDataHandle<GFS,Z> adddh(gfs,z);
  	if (gfs.gridView().comm().size()>1)
    	gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);


  	// Assemble constraints
  	typedef typename GFS::template
    	ConstraintsContainer<RF>::Type CC;
  	CC cc;
  	//Dune::PDELab::constraints(b,gfs,cc); // assemble constraints
  	
	// Make a local operator
  	typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;
  	LOP lop(problem);

  	// Make a global operator
  	typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  	//int degree = ptree.get("fem.degree",(int)1);
  	MBE mbe((int)pow(1+2*degree,dim));
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GO;
  	GO go(gfs,cc,gfs,cc,lop,mbe);



  	typedef massFEM<Problem<RF>,FEM> LOPm;
  	LOPm lopm(problem);

  	// Make a global operator
  	//typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  	//int degree = ptree.get("fem.degree",(int)1);
  	//MBE mbe((int)pow(1+2*degree,dim));
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOPm,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GOm;
  	GOm gom(gfs,cc,gfs,cc,lopm,mbe);


  	// Make instationary grid operator
  	/*typedef NonlinearHeatFEM<Problem<RF>,FEM> LOP;
  	LOP lop(problem);
  	typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  	//int degree = ptree.get("fem.degree",(int)1);
  	MBE mbe((int)pow(1+2*degree,dim));
  	typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,
                                     RF,RF,RF,CC,CC> GO0;
  	GO0 go0(gfs,cc,gfs,cc,lop,mbe);*/



  	// make coefficent Vectors
  	using X = Dune::PDELab::Backend::Vector<GFS,double>;
  	X x(gfs,0.0);

  	// represent operator as a matrix
  	//typedef typename Dune::PDELab::Backend::Matrix<GFS,GFS,RF, RF> M;//GO::Jacobian M;
	typedef typename GO::template MatrixContainer<RF>::Type M;
  	M m(go);
  	//std::cout << m.patternStatistics() << std::endl;
  	m = 0.0;
  	go.jacobian(x,m);

	Dune::PDELab::Backend::native(m)[0][0][0][0] = 1;
	
	Dune::PDELab::Backend::native(m)[0][1][0][0] = 0;


	Dune::PDELab::Backend::native(m)[Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1] = 1;
	Dune::PDELab::Backend::native(m)[Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-2][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1] = 0;

	Dune::PDELab::Backend::native(m)[0][0];
	for(int i=0; i<Dune::PDELab::Backend::native(m).M(); i++){
		for(int j=0; j<Dune::PDELab::Backend::native(m).N(); j++){ 
			if (Dune::PDELab::Backend::native(m).exists(i,j)) {
				std::cout << Dune::PDELab::Backend::native(m)[i][j][0][0]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl;
	}
	//std::ofstream matrix("Matrix");
	//Dune::printmatrix(matrix, m.base(), "M", "r", 6, 3);

	
    	//Dune::printmatrix(std::cout, m, "A", "row");
	//for(auto r =m.begin(); r !=m.end(); ++r){std::cout << *r <<std::endl;}

	//std::ostringstream name;
    	//name<<0<<": row";

	//m.base()[0][0][0][0];
    	//Dune::printmatrix(std::cout, m, "A", "row");

	//Dune::PDELab::Backend::native(m).mv(initial, z);



  	// Select a linear solver backend NEW IN PARALLEL
  	typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;
  	int verbose=0;
  	if (gfs.gridView().comm().rank()==0) verbose=1;
  	LS ls(gfs,100,verbose);


	//typedef Dune::PDELab::OverlappingOperator<GFS, decltype(Dune::PDELab::Backend::native(m)), X, X> OP;
	//OP ovlpOperator(gfs, Dune::PDELab::Backend::native(m)); 

	double dt = 0.02;
	for (int i=0; i<10; i++){
		ls.apply(m, sol, z, 0);
		for(auto r =sol.begin(); r !=sol.end(); ++r){std::cout << "vector " << *r <<std::endl;}
		//M.mv(sol, z);		
		//go.apply(sol, z);
		gom.jacobian_apply(sol, z);
		for(auto r =z.begin(); r !=z.end(); ++r){std::cout << "new rhs " << *r <<std::endl;}
		Dune::PDELab::Backend::native(z)[0][0]=0;
		Dune::PDELab::Backend::native(z)[Dune::PDELab::Backend::native(m).M()-1][0]=0;
		//std::exit(0);
	}

	// M um+1 = M um + dt A um+1
	// (M + dt A) um+1 = M um 




	t2 = MPI_Wtime(); 
	printf( "Elapsed time is %f\n", t2 - t1 ); 

	for(int i=0; i<Dune::PDELab::Backend::native(m).N(); i++)
				std::cout << "initial " << Dune::PDELab::Backend::native(initial)[i][0] << " num " << Dune::PDELab::Backend::native(sol)[i][0] << " alg " << Dune::PDELab::Backend::native(vgl)[i][0] << std::endl;

	//for(auto r =sol.begin(); r !=sol.end(); ++r){std::cout << "vector" << *r <<std::endl;}
  	// solve nonlinear problem
  	/*Dune::PDELab::Newton<GO,LS,Z> newton(go,z,ls);
  	newton.setReassembleThreshold(0.0);
  	newton.setVerbosityLevel(2);
  	newton.setReduction(1e-10);
  	newton.setMinLinearReduction(1e-4);
  	newton.setMaxIterations(25);
  	newton.setLineSearchMaxIterations(10);
  	newton.apply();*/



  	// Write VTK output file
  	Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
  	typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  	vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
                                         VTKF(zdgf,"fesol")));
  	vtkwriter.write("output",
        Dune::VTK::appendedraw);



  // Make a local operator
  //typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;
  //LOP lop(problem);

  // Make a global operator
  //typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  //int degree = ptree.get("fem.degree",(int)1);
  //MBE mbe((int)pow(1+2*degree,dim));
  //typedef Dune::PDELab::GridOperator<
  //  GFS,GFS,  /* ansatz and test space */
  //  LOP,      /* local operator */
  //  MBE,      /* matrix backend */
  //  RF,RF,RF, /* domain, range, jacobian field type*/
  //  CC,CC     /* constraints for ansatz and test space */
  //  > GO;
  //GO go(gfs,cc,gfs,cc,lop,mbe);





	return 0;
}
