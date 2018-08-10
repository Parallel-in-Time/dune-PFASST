#ifndef _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_


#include <memory>
#include <type_traits>

using std::shared_ptr;
using std::vector;

#include <vector>

#include <pfasst/sweeper/FE_impl.hpp>
#include <pfasst/contrib/fft.hpp>

#include "fe_manager.hpp"

//create mass + dt stiffness 
#include"nonlinearpoissonfem.hh"
#include"Lnonlinearpoissonfem.hh"
//#ifndef PI
//#define PI 3.1415926535897932385
//#endif
template<typename Number>
class Problem
{
public:
  typedef Number value_type;
  const double dt=0;
  const double t=0;	
  Problem () {}
  Problem<Number> (double dt_, double t_ =0 ): dt(dt_), t(t_) {}
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=1;
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
        return solution * std::exp(-t * dim * PI*PI*0.1);
  }
};

template<typename Number>
class LProblem
{
public:
  typedef Number value_type;
  const double nu=0;	
  LProblem () {}
  LProblem<Number> (double nu_): nu(nu_) {}
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=1;
	const double t=0;	
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
        return solution * std::exp(-t * dim * PI*PI*0.1);
  }
};

namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
        /*template<
                class EncapsulationTraits,
                size_t base_order,
                size_t dim,
                class... Ts
        >
        struct dune_sweeper_traits : public sweeper_traits<EncapsulationTraits,  Ts...>
        {

            using encap_traits = EncapsulationTraits;
            using encap_t = encap::Encapsulation<EncapsulationTraits>;
            using time_t = typename encap_traits::time_t;
            using spatial_t = typename encap_traits::spatial_t;

            //static constexpr size_t BASE_ORDER = base_order;
            //static constexpr size_t DIM = dim;
        };*/

      template<
        class SweeperTrait,
	class Mass,
        typename Enabled = void
      >
      class Heat_FE
        : public IMEX<SweeperTrait, Mass, Enabled>
      {
        /*static_assert(std::is_same<
                        typename SweeperTrait::encap_t::traits::dim_t,
                        std::integral_constant<size_t, 1>
                      >::value,
                      "Heat_FE Sweeper requires 2D data structures");*/

        public:
          using traits = SweeperTrait;

          static void init_opts();

          int finer;

        private:
          using spatial_t = typename traits::spatial_t;

          typename traits::time_t                        _t0{0.0};
          spatial_t                                      _nu{0.1};//0.1
          pfasst::contrib::FFT<typename traits::encap_t> _fft;
          vector<vector<spatial_t>>                      _lap;


          std::shared_ptr<fe_manager> FinEl;
	  

	  typedef typename SweeperTrait::encap_traits::mass_t GOO;
	  shared_ptr<typename SweeperTrait::encap_traits::mass_t> M_dtA_dune;
          shared_ptr<typename SweeperTrait::encap_traits::mass_t> test_dune;

	  typedef typename GOO::template MatrixContainer<double>::Type MM;
	  shared_ptr<MM> mm;
          //MatrixType M_dune;
          //MatrixType A_dune;


 	const static int dim=1;
	const int degree =1;
	int nelements;
	int nlevel;

        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        //Dune::FieldVector<double,dim> h = {1};	      
	//std::array<int,dim> n;
	//std::fill(n.begin(), n.end(), nelements);
        std::shared_ptr<Grid> gridp;// = std::shared_ptr<Grid>(new Grid(h,n));

        //gridp->refineOptions(false); // keep overlap in cells
        //typedef Grid::LeafGridView GV;
	typedef Grid::LevelGridView GV;
        //GV gv=gridp->leafGridView();
	//typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
	typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF,double,  1> FEM;  
	//typedef	Dune::PDELab::PkQkLocalFiniteElementMap<GV, DF, double, 1 > FEM;        


	 std::shared_ptr<FEM> fem; 
	//FEM fem(gridp->leafGridView()); //anststatt gv

  	// Make grid function space
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	 std::shared_ptr<GFS> gfs; //(gridp->leafGridView(),fem);


	typedef double RF; 
  	std::shared_ptr<Problem<RF>> problem; //(0.05);
  	std::shared_ptr<Problem<RF>> mass_problem;//(0.0);
  	std::shared_ptr<LProblem<RF>> laplace_problem;//(0.0);
  	// Assemble constraints
  	typedef typename GFS::template
    	ConstraintsContainer<RF>::Type CC;
  	CC cc1;
  	CC cc2;
	CC cc3;


	// Make a local operator
  	typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;
  	typedef LNonlinearPoissonFEM<LProblem<RF>,FEM> LLOP;
  	std::shared_ptr<LOP> lop; //(problem);
	std::shared_ptr<LOP> mass_lop;
	std::shared_ptr<LLOP> laplace_lop;
  	// Make a global operator
  	typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  	std::shared_ptr<MBE> mbe;
	//MBE mbe((int)pow(1+2*degree,dim));
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GO;
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LLOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> LGO;

	shared_ptr<LGO> A_dune;	

	//this->M_dune = std::make_shared<GO>(gfs,cc,gfs,cc,lop,mbe)


  	//LOP lopm(mass_problem);

	//this->M_dtA_dune = std::make_shared<GO>(gfs,cc,gfs,cc,lopm,mbe);
	//this->test_dune = std::make_shared<GO>(gfs,cc,gfs,cc,lop,mbe);


  	using X = Dune::PDELab::Backend::Vector<GFS,double>;
  	shared_ptr<X> x;//(gfs,0.0);
  	shared_ptr<X> y;//(gfs,0.0);

  	//represent operator as a matrix
	typedef typename GO::template MatrixContainer<RF>::Type M;
	std::shared_ptr<M> m;
	std::shared_ptr<M> mM;
  	//M m(*(this->M_dune));
	//mm = std::make_shared<M>(*(this->M_dune));

  	//m = 0.0;
  	//this->M_dune->jacobian(x,m);


  	// Select a linear solver backend NEW IN PARALLEL
  	typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;
  	//int verbose=0;
  	shared_ptr<Dune::PDELab::ISTLBackend_CG_AMG_SSOR<typename SweeperTrait::encap_traits::mass_t>> ls;
  	//if (gfs.gridView().comm().rank()==0) verbose=1;






          std::shared_ptr<GridType> grid;



          std::shared_ptr<BasisFunction> basis;



        //protected:
      public:  
          /*virtual shared_ptr<typename SweeperTrait::encap_t>
          evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;*/

          virtual shared_ptr<typename SweeperTrait::encap_t> 
          evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;

          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                      shared_ptr<typename SweeperTrait::encap_t> u,
                                      const typename SweeperTrait::time_t& t,
                                      const typename SweeperTrait::time_t& dt,
                                      const shared_ptr<typename SweeperTrait::encap_t> rhs) override;

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_error(const typename SweeperTrait::time_t& t);

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                 const typename SweeperTrait::time_t& t);

        public:


          explicit Heat_FE(std::shared_ptr<fe_manager>, size_t);


          Heat_FE(const Heat_FE<SweeperTrait, Mass, Enabled>& other)= default;
          Heat_FE(Heat_FE<SweeperTrait, Mass, Enabled>&& other)= default;
          virtual ~Heat_FE() = default;
          Heat_FE<SweeperTrait, Mass, Enabled>& operator=(const Heat_FE<SweeperTrait, Mass, Enabled>& other) = default;
          Heat_FE<SweeperTrait, Mass, Enabled>& operator=(Heat_FE<SweeperTrait, Mass, Enabled>&& other) = default;

          virtual void set_options() override;

          template<typename Basis>
          void assemble(Basis &basis);

          void assemble();

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);

          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;

          shared_ptr<GridType> get_grid() const;
      };
    }  // ::pfasst::examples::heat_FE
  }  // ::pfasst::examples
}  // ::pfasst

#include "FE_sweeper_impl.hpp"

#endif // _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
