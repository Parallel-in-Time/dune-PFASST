#include "FE_sweeper.hpp"

//#include "assemble.hpp"

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




        auto exact_solution = [](const auto  &x){
            double solution=1.0;
            for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
            return solution * std::exp(0);
        };


      template<class SweeperTrait, class MassType, typename Enabled>
      void
      Heat_FE<SweeperTrait, MassType, Enabled>::init_opts()
      {
        /*config::options::add_option<size_t>("Heat FE", "num_dofs",
                                            "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat FE", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat FE", "nu",
                                               "thermal diffusivity");*/
      }




        template<class SweeperTrait, class MassType, typename Enabled>
      Heat_FE<SweeperTrait, MassType, Enabled>::Heat_FE(size_t nlevel) // std::shared_ptr<fe_manager> FinEl,
        :   IMEX<SweeperTrait, MassType, Enabled>()

      {

	//setup the grid
	int nelements=10;

        Dune::FieldVector<double,dim> h = {1};	      
	std::array<int,dim> n;
	std::fill(n.begin(), n.end(), nelements);
        std::shared_ptr<GridType> gridp = std::shared_ptr<GridType>(new GridType(h,n));

        gridp->refineOptions(false); // keep overlap in cells
        GV gv=gridp->leafGridView();
        FEM fem(gv);

  	//GFS gfs(gv,fem);
	gfs = std::make_shared<GFS>(gv, fem);


  	/*using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  	Z z(gfs); // mass times initial value
	Z initial(gfs); // initial value
	Z vgl(gfs); // analytic solution at the end point
	Z sol(gfs);*/ // numeric solution

  	// Make a grid function out of it
  	//typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  	//ZDGF zdgf(gfs,z);

  	/*Problem<RF> problem;
  	auto glambda = [&](const auto& x){return problem.g(x);};
  	auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);

	Problem2<RF> problem2;
  	auto glambda2 = [&](const auto& x){return problem2.g(x);};
  	auto g2 = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda2); */

  	// Fill the coefficient vector
  	//Dune::PDELab::interpolate(g,gfs,initial);//z
  	//Dune::PDELab::interpolate(g2,gfs,vgl); // vgl = analytic solution at the end point 

	//initial = z;

	//for(auto r =z.begin(); r !=z.end(); ++r){std::cout << "vector" << *r <<std::endl;}
 
	

  	// make vector consistent NEW IN PARALLEL
  	/*Dune::PDELab::istl::ParallelHelper<GFS> grid_helper(gfs);
  	grid_helper.maskForeignDOFs(z);
  	Dune::PDELab::AddDataHandle<GFS,Z> adddh(gfs,z);
  	if (gfs.gridView().comm().size()>1){
    		gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
	}*/


  	// Assemble constraints
  	typedef typename GFS::template
    	ConstraintsContainer<RF>::Type CC;
  	CC cc;
  	
	// Make a local operator
  	typedef NonlinearPoissonFEM<FEM> LOP;
  	LOP lop(0);

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
  	GO go(*gfs,cc,*gfs,cc,lop,mbe);



  	typedef massFEM<FEM> LOPm;
  	LOPm lopm(0);

  	// Make a global operator
  	typedef Dune::PDELab::GridOperator<
    		GFS,GFS,  // ansatz and test space 
    		LOPm,      // local operator 
    		MBE,      // matrix backend 
    		RF,RF,RF, // domain, range, jacobian field type
    		CC,CC     // constraints for ansatz and test space 
    		> GOm;
  	GOm gom(*gfs,cc,*gfs,cc,lopm,mbe);

	//gom.jacobian_apply(initial, z);

  	// make coefficent Vectors
  	using X = Dune::PDELab::Backend::Vector<GFS,double>;
  	X x(*gfs,0.0);

  	// represent operator as a matrix
	typedef typename GO::template MatrixContainer<RF>::Type M;
  	M m(go);
  	//std::cout << m.patternStatistics() << std::endl;
  	m = 0.0;
  	go.jacobian(x,m);

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
  	typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;
  	int verbose=0;
  	if (gfs->gridView().comm().rank()==0) verbose=1;
  	LS ls(*gfs,100,verbose);




	//this->FinEl = FinEl;
	//basis = FinEl->get_basis(nlevel);
	    
	//assembleProblem(basis, this->A_dune, this->M_dune);

        const auto bs = Dune::PDELab::Backend::native(m).M(); //basis->size();
        std::cout << "Finite Element basis consists of " <<  bs << " elements " << std::endl;

        this->encap_factory()->set_size(bs);
        this->encap_factory()->set_gfs(*gfs);

      }




      template<class SweeperTrait, class MassType, typename Enabled>
      void
      Heat_FE<SweeperTrait, MassType, Enabled>::set_options()
      {


        IMEX<SweeperTrait, MassType, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        std::cout << "thermal diffusivity is set to "   << this->_nu << std::endl;

        int num_nodes = this->get_quadrature()->get_num_nodes();



      }

        /*template<class SweeperTrait, class MassType, typename Enabled>
        template<typename Basis>
      void
      Heat_FE<SweeperTrait, MassType, Enabled>::assemble(Basis &basis){

        assembleProblem(basis, this->A_dune, this->M_dune);


      };

      template<class SweeperTrait, class MassType, typename Enabled>
      void
      Heat_FE<SweeperTrait, MassType, Enabled>::assemble(){

        assembleProblem(basis, this->A_dune, this->M_dune);


      };*/




      template<class SweeperTrait, class MassType, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, MassType, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        //const auto dim = dim;
        spatial_t nu = this-> _nu;


        //interpolate(*gfs, result->data(), exact_solution);

	//std::exit(0);

        return result;
      }

      template<class SweeperTrait, class MassType, typename Enabled>
      void
      Heat_FE<SweeperTrait, MassType, Enabled>::post_step()
      {
        IMEX<SweeperTrait, MassType, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        //ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        //this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, class MassType, typename Enabled>
      bool
      Heat_FE<SweeperTrait, MassType, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, MassType, Enabled>::converged(pre_check);

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

      template<class SweeperTrait, class MassType, typename Enabled>
      bool
      Heat_FE<SweeperTrait, MassType, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, class MassType, typename Enabled>
      size_t
      Heat_FE<SweeperTrait, MassType, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }



      /*template<class SweeperTrait, class MassType, typename Enabled>
      shared_ptr<GridType>
      Heat_FE<SweeperTrait, MassType, Enabled>::get_grid() const
      {
        return grid;
      }*/


      template<class SweeperTrait, class MassType, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, MassType, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
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

      template<class SweeperTrait, class MassType, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, MassType, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
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

      /*template<class SweeperTrait, class MassType, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, MassType, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(),  "evaluating EXPLICIT part at t=" << t);
        auto result = this->get_encap_factory().create();
        result->zero();
        this->_num_expl_f_evals++;
        return result;
      }*/

      template<class SweeperTrait, class MassType, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, MassType, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {



        ML_CVLOG(4, this->get_logger_id(), "evaluating IMPLICIT part at t=" << t);


        auto result = this->get_encap_factory().create();

        /*double nu =this->_nu;

        this->A_dune.mmv(u->get_data(), result->data());


        result->data() *= nu;

        for (size_t i = 0; i < u->get_data().size(); i++) {

        }*/

        
        return result;
        


      }

      template<class SweeperTrait, class MassType, typename Enabled>
      void
      Heat_FE<SweeperTrait, MassType, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {


        
        /*Dune::BlockVector<Dune::FieldVector<double,1> > M_rhs_dune ;
        M_rhs_dune.resize(rhs->get_data().size());      
	M_rhs_dune = rhs->get_data(); 
        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dtA_dune = 	Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->A_dune);
        M_dtA_dune *= (dt * this->_nu); // fehler ruth 
        M_dtA_dune += this->M_dune;
        Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dtA_dune);
        Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dtA_dune,1.0);
        Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-10, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver
        Dune::InverseOperatorResult statistics ;
        cg.apply(u->data(), M_rhs_dune , statistics ); //rhs ist nicht constant!!!!!!!!!
        Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
        this->M_dune.mv(u->get_data(), M_u);
        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
        }*/
        this->_num_impl_solves++;






      }
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
} // ::pfasst
