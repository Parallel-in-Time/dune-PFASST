#include "FE_sweeper.hpp"

#include "assemble.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <utility>
#include <vector>
using std::shared_ptr;
using std::vector;


#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>

#include <iostream>





namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {


/*template<typename Number>
class Problem
{
public:
  typedef Number value_type;
  Problem () {}
  template<typename X>
  Number g (const X& x) const
  {
	const int dim=1;
	const double t=0;	
	double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
        return solution * std::exp(-t * dim * PI*PI);
  }
};*/





      template<class SweeperTrait, class Mass, typename Enabled>
      void
      Heat_FE<SweeperTrait, Mass, Enabled>::init_opts()
      {
        /*config::options::add_option<size_t>("Heat FE", "num_dofs",
                                            "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat FE", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat FE", "nu",
                                               "thermal diffusivity");*/
      }


        template<class SweeperTrait, class Mass, typename Enabled>
      Heat_FE<SweeperTrait, Mass, Enabled>::Heat_FE( std::shared_ptr<fe_manager> FinEl, size_t nlevel) // std::shared_ptr<fe_manager> FinEl,
        :   IMEX<SweeperTrait, Mass, Enabled>()

      {


	this->nelements =FinEl->get_nelem();
        //this->encap_factory()->set_size(nelements);
	
        Dune::FieldVector<double,dim> h = {1};	      
	std::array<int,dim> n;
	std::fill(n.begin(), n.end(), nelements);
        gridp = std::shared_ptr<Grid>(new Grid(h,n));

        gridp->refineOptions(false); 


        fem = std::make_shared<FEM>(gridp->leafGridView());

  	// Make grid function space
  	gfs  = std::make_shared<GFS>(gridp->leafGridView(),*fem);


        this->encap_factory()->set_gfs(*gfs);
  	
	/*using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
	Z z(*gfs);
  	typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  	ZDGF zdgf(*gfs,z);
  	Dune::PDELab::istl::ParallelHelper<GFS> grid_helper(*gfs);
  	grid_helper.maskForeignDOFs(z);
  	Dune::PDELab::AddDataHandle<GFS,Z> adddh(*gfs,z);
  	if (gfs->gridView().comm().size()>1){
    		gfs->gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
	}*/


	typedef double RF; 
  	problem = std::make_shared<Problem<RF>>(1.0); //this->get_status()->get_dt() 0.00984077
  	mass_problem = std::make_shared<Problem<RF>>(0.0);
  	laplace_problem = std::make_shared<LProblem<RF>>(this->_nu);
  	// Assemble constraints
  	//typedef typename GFS::template
    	//ConstraintsContainer<RF>::Type CC;
  	//CC cc;
  	
	// Make a local operator
  	//typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;

  	lop = std::make_shared<LOP>(*problem);
  	mass_lop = std::make_shared<LOP>(*mass_problem);
  	laplace_lop = std::make_shared<LLOP>(*laplace_problem);

  	// Make a global operator
  	/*typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  	MBE mbe((int)pow(1+2*degree,dim));*/
	mbe = std::make_shared<MBE>((int)pow(1+2*degree,dim));
  	/*typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOP,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GO;*/

	this->M_dune = std::make_shared<GO>(*gfs,cc1,*gfs,cc1,*mass_lop,*mbe);
	this->M_dtA_dune = std::make_shared<GO>(*gfs,cc2,*gfs,cc2,*lop,*mbe);
	this->A_dune = std::make_shared<LGO>(*gfs,cc3,*gfs,cc3,*laplace_lop,*mbe);

	/*M m2(*(this->M_dune));
	X x2(*gfs, 0.0);
   	this->M_dune->jacobian(x2,m2);


	for(int i=0; i<Dune::PDELab::Backend::native(m2).M(); i++){
		for(int j=0; j<Dune::PDELab::Backend::native(m2).N(); j++){ 
			if (Dune::PDELab::Backend::native(m2).exists(i,j)) {
				std::cout << Dune::PDELab::Backend::native(m2)[i][j][0][0]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl;
	}std::exit(0);*/


	/*auto t1 = this->get_encap_factory().create();
	auto t2 = this->get_encap_factory().create();
	this->M_dune->jacobian_apply(t1->data(), t2->data());

	std::cout << "*******************************   multipliziert *********************** " <<  std::endl;*/
	//this->M_dtA_dune = std::make_shared<GO>(*gfs,cc,*gfs,cc,lopm,mbe);

  	/*LOP lopm(mass_problem);

	this->M_dtA_dune = std::make_shared<GO>(*gfs,cc,*gfs,cc,lopm,mbe);
	this->test_dune = std::make_shared<GO>(*gfs,cc,*gfs,cc,lop,mbe);
  	//GO gom(gfs,cc,gfs,cc,lopm,mbe);


  	using X = Dune::PDELab::Backend::Vector<GFS,double>;*/
  	this->x = std::make_shared<X>(*gfs,0.0);

  	//represent operator as a matrix
  	m = std::make_shared<M>(*(this->M_dtA_dune));
	//mm = std::make_shared<M>(*(this->M_dune));

  	*m = 0.0;
  	this->M_dtA_dune->jacobian(*x,*m);

	Dune::PDELab::Backend::native(*m)[0][0][0][0] = 1;
	
	Dune::PDELab::Backend::native(*m)[0][1][0][0] = 0;


	Dune::PDELab::Backend::native(*m)[Dune::PDELab::Backend::native(*m).M()-1][Dune::PDELab::Backend::native(*m).M()-1][Dune::PDELab::Backend::native(*m).M()-1][Dune::PDELab::Backend::native(*m).M()-1] = 1;
	Dune::PDELab::Backend::native(*m)[Dune::PDELab::Backend::native(*m).M()-1][Dune::PDELab::Backend::native(*m).M()-2][Dune::PDELab::Backend::native(*m).M()-1][Dune::PDELab::Backend::native(*m).M()-1] = 0;



	for(int i=0; i<Dune::PDELab::Backend::native(*m).M(); i++){
		for(int j=0; j<Dune::PDELab::Backend::native(*m).N(); j++){ 
			if (Dune::PDELab::Backend::native(*m).exists(i,j)) {
				std::cout << Dune::PDELab::Backend::native(*m)[i][j][0][0]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl;
	}//std::exit(0);
	//std::cout << "****************************************************** " << Dune::PDELab::Backend::native(*m)[1][1][1][1] << std::endl;

  	// Select a linear solver backend NEW IN PARALLEL
  	//typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;
  	int verbose=0;
  	if (gfs->gridView().comm().rank()==0) verbose=1;
  	this->ls = std::make_shared<LS>(*gfs,nelements,verbose);




      }




      template<class SweeperTrait, class Mass, typename Enabled>
      void
      Heat_FE<SweeperTrait, Mass, Enabled>::set_options()
      {


        IMEX<SweeperTrait, Mass, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        std::cout << "thermal diffusivity is set to "   << this->_nu << std::endl;

        int num_nodes = this->get_quadrature()->get_num_nodes();



      }

        template<class SweeperTrait, class Mass, typename Enabled>
        template<typename Basis>
      void
      Heat_FE<SweeperTrait, Mass, Enabled>::assemble(Basis &basis){

        assembleProblem(basis, this->A_dune, this->M_dune);


      };

      template<class SweeperTrait, class Mass, typename Enabled>
      void
      Heat_FE<SweeperTrait, Mass, Enabled>::assemble(){

        assembleProblem(basis, this->A_dune, this->M_dune);


      };




      template<class SweeperTrait, class Mass, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Mass, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {

	std::cout << "*******************************  im exact *********************** " <<  std::endl;

	auto t1 = this->get_encap_factory().create();
	auto t2 = this->get_encap_factory().create();
	this->M_dune->jacobian_apply(t1->data(), t2->data());



        auto result = this->get_encap_factory().create();

        const auto dim = DIMENSION;
        spatial_t nu = this-> _nu;


	const int degree =1;
        typedef Dune::YaspGrid<1> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<double,1> h = {1};	      
	std::array<int,1> n;
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

  	Problem<double> problem(0,t);
  	auto glambda = [&](const auto& x){return problem.g(x);};
  	auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);

 	using Z = Dune::PDELab::Backend::Vector<GFS,double>;
  	Z z(gfs); // mass times initial value
	//std::cout << "vor " << std::endl;
	//std::cout << "vor " << result->get_total_num_dofs() << std::endl;
  	Dune::PDELab::interpolate(g,gfs, result->data());//z
	//std::cout << "nach " << z.N() << std::endl;
	//std::cout << "nach " << std::endl;
	//result->data() = z;
        /*auto exact_solution = [t, nu, dim](const InVectorType &x){
            double solution=1.0;
            for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
            return solution * std::exp(-t * dim * PI*PI * nu);
        };

        

        //interpolate(*basis, result->data(), exact_solution);
       for (size_t i = 0; i < result->get_data().N(); i++) {
          //std::cout << "result exact" << result->data()[i] << std::endl;
        }
	//std::exit(0);*/



        return result; //std::make_shared<typedef(z)>(z);
      }

      template<class SweeperTrait, class Mass, typename Enabled>
      void
      Heat_FE<SweeperTrait, Mass, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Mass, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        //ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        //this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, class Mass, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Mass, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, Mass, Enabled>::converged(pre_check);

        if (!pre_check) {
          assert(this->get_status() != nullptr);
          const typename traits::time_t t = this->get_status()->get_time();
          const typename traits::time_t dt = this->get_status()->get_dt();

          assert(this->get_quadrature() != nullptr);
          auto nodes = this->get_quadrature()->get_nodes();
          const auto num_nodes = this->get_quadrature()->get_num_nodes();
          nodes.insert(nodes.begin(), typename traits::time_t(0.0));

          ML_CVLOG(1, this->get_logger_id(),
                   "Observables after "
                   << ((this->get_status()->get_iteration() == 0)
                          ? std::string("prediction")
                          : std::string("iteration ") + std::to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(),
                     "  t["<<m<<"]=" <<  (t + dt * nodes[m])
                     << "      |abs residual| = " <<  this->_abs_res_norms[m]
                     << "      |rel residual| = " <<  this->_rel_res_norms[m]
//                      << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
//                      << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" <<  (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " << this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " << this->_rel_res_norms[num_nodes]
//                   << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
//                   << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes])
                 );
        }
        return converged;
      }

      template<class SweeperTrait, class Mass, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Mass, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, class Mass, typename Enabled>
      size_t
      Heat_FE<SweeperTrait, Mass, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().N();
      }



      template<class SweeperTrait, class Mass, typename Enabled>
      shared_ptr<GridType>
      Heat_FE<SweeperTrait, Mass, Enabled>::get_grid() const
      {
        return grid;
      }


      template<class SweeperTrait, class Mass, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Mass, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
      {
        ML_CVLOG(4, this->get_logger_id(), "computing error");

        assert(this->get_status() != nullptr);
        const typename traits::time_t dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), typename traits::time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> error;
        error.resize(num_nodes + 1);
        std::generate(error.begin(), error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const typename traits::time_t ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), this->get_states()[m]);
        }

        return error;
      }

      template<class SweeperTrait, class Mass, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Mass, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                            const typename SweeperTrait::time_t& t)
      {
        UNUSED(t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> rel_error;
        rel_error.resize(error.size());
        std::generate(rel_error.begin(), rel_error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          rel_error[m]->scaled_add(1.0 / this->get_states()[m]->norm0(), error[m]);
        }

        return rel_error;
      }



      template<class SweeperTrait, class Mass, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Mass, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 start evaluate" <<  std::endl;
	//std::cout << "im evaluate ACHTUNG!!!   " << std::endl; std::exit(0);

        ML_CVLOG(4, this->get_logger_id(), "evaluating IMPLICIT part at t=" << t);


        auto f = this->get_encap_factory().create();
	this->A_dune->jacobian_apply((u->data()), (f->data()));

	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

	if(my_rank==num_pro-1) Dune::PDELab::Backend::native(f->data())[f->data().N()-1][0] = 0;
	if(my_rank==0) Dune::PDELab::Backend::native(f->data())[0][0] = 0;

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ende evaluate" <<  std::endl;
	return f;
        //f->data() = result->data(); 
        //return std::make_shared<typename SweeperTrait::encap_t>(*result);
        


      }

      /*template<class SweeperTrait, class Mass, typename Enabled>
      typename SweeperTrait::encap_t
      Heat_FE<SweeperTrait, Mass, Enabled>::evaluate_rhs_impl2(const typename SweeperTrait::time_t& t,
                                                       const typename SweeperTrait::encap_t u)
      {

	//std::cout << "im evaluate ACHTUNG!!!   " << std::endl; std::exit(0);

        ML_CVLOG(4, this->get_logger_id(), "evaluating IMPLICIT part at t=" << t);


        auto result = this->get_encap_factory().create();
	this->A_dune->jacobian_apply(u->data(), result->data());



	Dune::PDELab::Backend::native(result->data())[100][0] = 0;
	Dune::PDELab::Backend::native(result->data())[0][0] = 0;

	//for(auto r =result->data().begin(); r !=result->data().end(); ++r){std::cout << "evaluate " << *r <<std::endl;} std::exit(0);
        
        return (*result);
        


      }*/

      template<class SweeperTrait, class Mass, typename Enabled>
      void
      Heat_FE<SweeperTrait, Mass, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {


	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 start solve" << std::endl;
	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

	if (my_rank == 0) for(auto r =u->data().begin(); r !=u->data().end(); ++r){std::cout << my_rank << " " << u->data().N() << " evaluate " << *r <<std::endl;} 
        MPI_Barrier(MPI_COMM_WORLD);
	if (my_rank == 1) for(auto r =u->data().begin(); r !=u->data().end(); ++r){std::cout << my_rank << " " << u->data().N() << " evaluate " << *r <<std::endl;} //std::exit(0);



	auto M_rhs_dune = this->get_encap_factory().create();
	std::cout << "im implicit solve " << std::endl;  


	if (my_rank==num_pro-1) Dune::PDELab::Backend::native(rhs->data())[rhs->data().N()-1][0] = 0;
	if(my_rank==0) Dune::PDELab::Backend::native(rhs->data())[0][0] = 0;

	M_rhs_dune->data() = rhs->get_data();
	std::cout << "vorm loesen  " << std::endl;  

	//for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << "M result " << *r <<std::endl;} std::exit(0);

 	MPI_Barrier(MPI_COMM_WORLD);
	//std::cout << dt << std::endl; std::exit(0);
  	Problem<RF> mdta_problem(dt * this->_nu); 
  	

  	LOP mdta_lop(mdta_problem);
  	MBE mbe((int)pow(1+2*degree,dim));
	CC cc4;
	GO mdta_go(*gfs,cc4,*gfs,cc4,mdta_lop,mbe); 
	

  	X x(*gfs,0.0);

  	M m(mdta_go);


  	m = 0.0;
  	mdta_go.jacobian(x,m);

 	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "bevor ich die randwerte setze  " << std::endl;  
	if(my_rank==0 ) Dune::PDELab::Backend::native(m)[0][0][0][0] = 1;
	
	if(my_rank==0 ) Dune::PDELab::Backend::native(m)[0][1][0][0] = 0;


	if(my_rank==num_pro-1 ) Dune::PDELab::Backend::native(m)[Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1] = 1;
	if(my_rank==num_pro-1 ) Dune::PDELab::Backend::native(m)[Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-2][Dune::PDELab::Backend::native(m).M()-1][Dune::PDELab::Backend::native(m).M()-1] = 0;
	std::cout << "nachdem ich die randwerte setze  " << std::endl;  


 	MPI_Barrier(MPI_COMM_WORLD);


	std::cout << "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  " <<my_rank << std::endl;  
	if(my_rank==0 )for(int i=0; i<Dune::PDELab::Backend::native(m).M(); i++){
		for(int j=0; j<Dune::PDELab::Backend::native(m).N(); j++){ 
			if (Dune::PDELab::Backend::native(m).exists(i,j)) {
				std::cout << Dune::PDELab::Backend::native(m)[i][j][0][0]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl;
	}        MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  " <<my_rank << std::endl;  
	if(my_rank==1 )for(int i=0; i<Dune::PDELab::Backend::native(m).M(); i++){
		for(int j=0; j<Dune::PDELab::Backend::native(m).N(); j++){ 
			if (Dune::PDELab::Backend::native(m).exists(i,j)) {
				std::cout << Dune::PDELab::Backend::native(m)[i][j][0][0]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl;
	}


	//std::cout << "bevor ich den quatsch hier lese" << std::endl;
	this->ls->apply(m, u->data(), M_rhs_dune->data(), 0);
	std::cout << "im impl solve " << std::endl;
	//for(auto r =u->data().begin(); r !=u->data().end(); ++r){std::cout << "Gleichungssystem geloest " << *r <<std::endl;}
	//std::exit(0);
        /*ML_CVLOG(4, this->get_logger_id(),
                 "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);


	
        Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().N());*/


	std::cout << "bevor ich den quatsch hier lese" << std::endl; 

	//for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << "f data " << *r <<std::endl;}std::exit(0);

	auto M_u = this->get_encap_factory().create();
        //*this->M_dune.mv(u->get_data(), M_u->data());
	//M_u->apply_Mass(this->M_dune, u);
	u->apply_Mass(this->M_dune, M_u);
	//M_dune->jacobian_apply(u->data(), );
	
	//for(auto r =M_u->data().begin(); r !=M_u->data().end(); ++r){std::cout << "M u data " << *r <<std::endl;}std::exit(0);

	M_u->data() -= rhs->get_data();
	M_u->data() *= 1./(dt);

	f->data() = M_u->data();



	if (my_rank==num_pro-1) Dune::PDELab::Backend::native(f->data())[f->data().N()-1][0] = 0;
	if (my_rank==0) Dune::PDELab::Backend::native(f->data())[0][0] = 0;


	//for(auto r =f->data().begin(); r !=f->data().end(); ++r){std::cout << "f data " << *r <<std::endl;}std::exit(0);
        this->_num_impl_solves++;

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ende solve" <<  std::endl;




      }
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
} // ::pfasst
