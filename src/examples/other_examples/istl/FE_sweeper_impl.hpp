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

	this->FinEl = FinEl;
	basis = FinEl->get_basis(nlevel);

        std::cout << "assemble " << std::endl;
	(this->A_dune)= std::make_shared<typename SweeperTrait::encap_traits::mass_t>();
        std::cout << "assemble " << std::endl;
	(this->M_dune)= std::make_shared<typename SweeperTrait::encap_traits::mass_t>();    
	assembleProblem(basis, *(this->A_dune), *(this->M_dune));

        const auto bs = basis->size();
        std::cout << "Finite Element basis consists of " <<  basis->size() << " elements " << std::endl;

        this->encap_factory()->set_size(bs);

	//setup the grid
        /*const int dim=1;
	const int degree =1;
	const int nelements =20;
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



        this->encap_factory()->set_gfs(gfs);*/


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
        auto result = this->get_encap_factory().create();

        const auto dim = DIMENSION;
        spatial_t nu = this-> _nu;

        auto exact_solution = [t, nu, dim](const InVectorType &x){
            double solution=1.0;
            for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
            return solution * std::exp(-t * dim * PI*PI * nu);
        };

        auto N_x = [t](const InVectorType &x){
            
            return x;

        };

        VectorType x_node;
        interpolate(*basis, x_node, N_x);

        interpolate(*basis, result->data(), exact_solution);
       for (size_t i = 0; i < result->get_data().size(); i++) {
          //std::cout << "result exact" << result->data()[i] << std::endl;
        }
	//std::exit(0);

        return result;
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
        return this->get_encap_factory().size();
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

      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(),  "evaluating EXPLICIT part at t=" << t);
        auto result = this->get_encap_factory().create();
        result->zero();
        this->_num_expl_f_evals++;
        return result;
      }*/

      template<class SweeperTrait, class Mass, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Mass, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 start evaluate" << std::endl;
	//std::cout << "im evaluate ACHTUNG!!!   " << std::endl; std::exit(0);
        ML_CVLOG(4, this->get_logger_id(), "evaluating IMPLICIT part at t=" << t);

        auto result = this->get_encap_factory().create();

        double nu =this->_nu;

        (*(this->A_dune)).mmv(u->get_data(), result->data());


        result->data() *= nu;
        //std::cout << "f_impl mit evaluate " << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          //f->data()[i] = (u->data()[i] - rhs->data()[i]) / (dt);
          //std::cout << "f u " << result->data()[i] << std::endl;
        }

        	//for(auto r =result->data().begin(); r !=result->data().end(); ++r){std::cout << "evaluate " << *r <<std::endl;} std::exit(0);
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 end evaluate" << std::endl;
        return result;
        


      }

      template<class SweeperTrait, class Mass, typename Enabled>
      void
      Heat_FE<SweeperTrait, Mass, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {

       std::cout << "das ist jetzt das rhs " <<  std::endl;
       for (size_t i = 0; i < rhs->get_data().size(); i++) {
          std::cout << "evaluate " << rhs->data()[i] << std::endl;
        }
	//std::exit(0);

	//if (my_rank == 0) for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << my_rank << " " << rhs->data().N() << " evaluate " << *r <<std::endl;} 
        //MPI_Barrier(MPI_COMM_WORLD);
	//if (my_rank == 1) for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << my_rank << " " << rhs->data().N() << " evaluate " << *r <<std::endl;} std::exit(0);

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 start impl solve" <<  std::endl;
        
        Dune::BlockVector<Dune::FieldVector<double,1> > M_rhs_dune ;
        M_rhs_dune.resize(rhs->get_data().size());
	
        
	M_rhs_dune = rhs->get_data(); 

	

	//std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 dt " << dt<< std::endl;
	//std::exit(0);
        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dtA_dune = 	Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(*this->A_dune);
	//dt=1;
	//std::cout <<  dt * this->_nu << " dt ";
        M_dtA_dune *= (dt * this->_nu); // fehler ruth 
        M_dtA_dune += *(this->M_dune);


	

        /*std::cout << "das ist jetzt die matrix dt " << dt <<  std::endl;	
	for (size_t i = 0; i < 101; i++) {
			for (size_t j = 0; j < 101; j++){
          			if (M_dtA_dune.exists(i,j)) {std::cout <<  ((M_dtA_dune))[i][j] << " ";}else{std::cout  << 0 << " ";} }std::cout << std::endl;
        } std::cout << dt * this->_nu << std::endl; std::exit(0);*/

        Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dtA_dune);

        Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dtA_dune,1.0);

        Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-10, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver

        Dune::InverseOperatorResult statistics ;

        cg.apply(u->data(), M_rhs_dune , statistics ); //rhs ist nicht constant!!!!!!!!!

          std::cout << "im impl solve  "  << std::endl;
       /*for (size_t i = 0; i < u->get_data().size(); i++) {
          std::cout << "Gleichungssystem geloest " << u->data()[i] << std::endl;
        }*/
	
	
	
	
        ML_CVLOG(4, this->get_logger_id(),
                 "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);


	
        Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
        (*(this->M_dune)).mv(u->get_data(), M_u);

        /*for (size_t i = 0; i < u->get_data().size(); i++) {
          //f->data()[i] = (u->data()[i] - rhs->data()[i]) / (dt);
          std::cout << "M u " << M_u[i] << std::endl;
        }std::exit(0);*/	

        //std::cout << "f_impl mit impl_solve" << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          //f->data()[i] = (u->data()[i] - rhs->data()[i]) / (dt);
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
          //std::cout << "M u " << f->data()[i] << std::endl;
        }//std::exit(0);

        this->_num_impl_solves++;
        //if (this->_num_impl_solves==5) std::exit(0);

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                 end impl solve" << std::endl;



      }
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
} // ::pfasst
