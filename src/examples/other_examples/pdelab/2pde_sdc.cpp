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
//create mass + dt stiffness 
#include"nonlinearpoissonfem.hh"
//creates mass
#include"mass_operator.hh"
//constants
#ifndef PI
#define PI 3.1415926535897932385
#endif
const int DIM=1;
const int nelements=100;


template<typename Number>
class Problem
{
public:
  typedef Number value_type;
  const double dt;	
  Problem () {}
  Problem<Number> (double dt_): dt(dt_) {}
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
  const double dt;
  typedef Number value_type;
  Problem2 () {}
  Problem2 (double dt_): dt(dt_) {}
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

	//setup the grid
        const int dim=DIM;
	const int degree =1;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
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

  	// Make grid function space
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	GFS gfs(gv,fem);

	typedef double RF; 
  	using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  	Z z(gfs); // mass times initial value
	Z initial(gfs); // initial value
	Z vgl(gfs); // analytic solution at the end point
	Z sol(gfs); // numeric solution

  	// Make a grid function out of it
  	typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  	ZDGF zdgf(gfs,z);

  	Problem<RF> problem(0.002);
  	auto glambda = [&](const auto& x){return problem.g(x);};
  	auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);

	Problem2<RF> problem2(0);
  	auto glambda2 = [&](const auto& x){return problem2.g(x);};
  	auto g2 = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda2); 

  	// Fill the coefficient vector
  	Dune::PDELab::interpolate(g,gfs,initial);//z
  	Dune::PDELab::interpolate(g2,gfs,vgl); // vgl = analytic solution at the end point 

	/*for(int i=0; i<nelements; i++)
				std::cout << "initial " << Dune::PDELab::Backend::native(initial)[i][0] << "   " << Dune::PDELab::Backend::native(vgl)[i][0] << std::endl;
	std::exit(0);*/

	//initial = z;

	//for(auto r =z.begin(); r !=z.end(); ++r){std::cout << "vector" << *r <<std::endl;}
 
	

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

  	typedef NonlinearPoissonFEM<Problem2<RF>,FEM> LOPm;
  	LOPm lopm(problem2);

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



  	//typedef massFEM<Problem<RF>,FEM> LOPm;
  	//LOP lopm(problem2);

  	// Make a global operator
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GOm;
  	GO gom(gfs,cc,gfs,cc,lopm,mbe);

	gom.jacobian_apply(initial, z);

  	// make coefficent Vectors
  	using X = Dune::PDELab::Backend::Vector<GFS,double>;
  	X x(gfs,0.0);

  	// represent operator as a matrix
	typedef typename GO::template MatrixContainer<RF>::Type M;
  	M mdta(go);
  	//std::cout << m.patternStatistics() << std::endl;
  	mdta = 0.0;
  	go.jacobian(x,mdta);

	Dune::PDELab::Backend::native(mdta)[0][0][0][0] = 1;
	
	Dune::PDELab::Backend::native(mdta)[0][1][0][0] = 0;


	Dune::PDELab::Backend::native(mdta)[Dune::PDELab::Backend::native(mdta).M()-1][Dune::PDELab::Backend::native(mdta).M()-1][Dune::PDELab::Backend::native(mdta).M()-1][Dune::PDELab::Backend::native(mdta).M()-1] = 1;
	Dune::PDELab::Backend::native(mdta)[Dune::PDELab::Backend::native(mdta).M()-1][Dune::PDELab::Backend::native(mdta).M()-2][Dune::PDELab::Backend::native(mdta).M()-1][Dune::PDELab::Backend::native(mdta).M()-1] = 0;

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
  	typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;
  	int verbose=0;
  	if (gfs.gridView().comm().rank()==0) verbose=1;
  	LS ls(gfs,100,verbose);


	for (int i=0; i<10; i++){
		ls.apply(mdta, sol, z, 0);
		//for(auto r =sol.begin(); r !=sol.end(); ++r){std::cout << "vector " << *r <<std::endl;}
		gom.jacobian_apply(sol, z);
		Dune::PDELab::Backend::native(z)[0][0]=0;
		Dune::PDELab::Backend::native(z)[Dune::PDELab::Backend::native(mdta).M()-1][0]=0;
		//for(auto r =z.begin(); r !=z.end(); ++r){std::cout << "new rhs " << *r <<std::endl;}
		//std::exit(0);
	}

	// M um+1 = M um + dt A um+1
	// (M + dt A) um+1 = M um 


	/*vgl -= sol; 
	for(auto r =vgl.begin(); r !=vgl.end(); ++r){std::cout << "vector " << *r <<std::endl;}
	//gom.applyscaleadd(-1, sol, vgl);
	std::cout << "norm    " << Dune::PDELab::Backend::native(vgl).infinity_norm() << std::endl;*/


	for(int i=0; i<Dune::PDELab::Backend::native(mdta).N(); i++)
				std::cout << "initial " << Dune::PDELab::Backend::native(initial)[i][0] << " num " << Dune::PDELab::Backend::native(sol)[i][0] << " alg " << Dune::PDELab::Backend::native(vgl)[i][0] << std::endl;

	t2 = MPI_Wtime(); 
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
