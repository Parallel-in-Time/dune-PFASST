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
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::init_opts()
      {

      }

      
      template<class SweeperTrait, typename Enabled>
      Heat_FE<SweeperTrait, Enabled>::Heat_FE(std::shared_ptr<fe_manager> FinEl, size_t nlevel)
        :   IMEX<SweeperTrait, Enabled>()

      {
      

	this->FinEl = FinEl;
	basis = FinEl->get_basis(nlevel);
	    
	assembleProblem(basis, this->A_dune, this->M_dune);

        const auto bs = basis->size();
        std::cout << "Finite Element basis of level " << nlevel << " consists of " <<  basis->size() << " elements " << std::endl;

        this->encap_factory()->set_size(bs);

      }
      

       
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::set_options()
      {

        IMEX<SweeperTrait, Enabled>::set_options();
        int num_nodes = this->get_quadrature()->get_num_nodes();


      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

	spatial_t n  = this-> _n;
    	spatial_t l0 = this-> _nu;
	spatial_t l1 = l0/2. *(pow((1+n/2.), 1/2.) + pow((1+ n/2.), -1/2.) );
	spatial_t d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);

        auto exact_solution = [l0, l1, n, d, t](const Dune::FieldVector<double,DIMENSION>&x){ 
	  return pow((1 + (pow(2, n/2.)-1 )* exp(-(n/2.)*d*(x+2*l1*t)) ), -2./n);
        };	

        auto N_x = [t](const Dune::FieldVector<double,DIMENSION>&x){
            return x;

        };

        Dune::BlockVector<Dune::FieldVector<double,DIMENSION>> x_node;
        interpolate(*basis, x_node, N_x);

        interpolate(*basis, result->data(), exact_solution);

        return result;
      }
      

      

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();
        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, Enabled>::converged(pre_check);

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
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" <<  (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " <<  this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " <<  this->_rel_res_norms[num_nodes]
                 );
        }
        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      Heat_FE<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }



      template<class SweeperTrait, typename Enabled> //Fehler der aktuellen Loesung an jedem Quadraturpunkt
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
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

      template<class SweeperTrait, typename Enabled> //vector encap_ raumdaten an jedem Quadraturpunkt
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
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



      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
	
         ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);



        auto result = this->get_encap_factory().create();
        auto u2 = this->get_encap_factory().create();
        double nu =this->_nu;


	u2->zero();
	for (int i=0; i<u->get_data().size(); ++i)
	    {
	    u2->data()[i]= -pow(u->get_data()[i], _n+1);	
	    }
	this->M_dune.mv(u2->get_data(), result->data());
	this->M_dune.umv(u->get_data(), result->data());
	result->data()*=(this->_nu*this->_nu);
	this->A_dune.umv(u->get_data(), result->data());

	

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {
	
        ML_CVLOG(4, this->get_logger_id(),
                 "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);

    	auto residuum = this->get_encap_factory().create();
	Dune::BlockVector<Dune::FieldVector<double,1> > newton_rhs;
    	newton_rhs.resize(rhs->get_data().size());
    
        int newton_max_steps=200;
	for (int i=0; i< newton_max_steps ;i++){
	  //compute jacobi-matrix and rhs in every Newton-step
	  Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune);
	  evaluate_f(f, u, dt, rhs);
	  evaluate_df(df, u, dt);
	  df.mv(u->data(), newton_rhs);
	  newton_rhs -= f->data();

	  Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df);
	  Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df,1);
	  Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-12, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver
	  Dune::InverseOperatorResult statistics ;
	  cg.apply(u->data(), newton_rhs , statistics ); //rhs ist nicht constant!!!!!!!!!
	  num_solves++;

          evaluate_f(f, u, dt, rhs);

          if(f->norm0()<this->newton){ if(!this->is_coarse) std::cout << "***************************************** anzahl iterationen innerer newton " << i+1 << " " << num_solves <<std::endl;   break;}
	  

	}	
	
	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);

        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
        }

        this->_num_impl_solves++;

      }
      
      
      template<class SweeperTrait, typename Enabled> void
      Heat_FE<SweeperTrait, Enabled>::evaluate_f(shared_ptr<typename SweeperTrait::encap_t> f,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt,
						 const shared_ptr<typename SweeperTrait::encap_t> rhs
						){
          
          
        f->zero();
	Dune::BlockVector<Dune::FieldVector<double,1> > fneu;
        fneu.resize(u->get_data().size());
	for (int i=0; i<u->get_data().size(); ++i)
	{
	  fneu[i]= (this->_nu*this->_nu) * (pow(u->get_data()[i], _n+1) - u->get_data()[i]);	
	}
	this->M_dune.mv(fneu, f->data());
	this->A_dune.mmv(u->get_data(),f->data());
	f->data() *= dt;
	this->M_dune.umv(u->get_data(),f->data());
	f->data() -=rhs->get_data();

      }

template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						){
            double _nu=this->_nu;
            double _n=this->_n;
            for (int i=0; i<df.N(); ++i)
            {
		for(int j=0; j< df.M(); j++)
                	if (df.exists(i,j)) df[i][j]= (_nu*_nu)*(_n+1) * pow(u->get_data()[j], _n) * this->M_dune[i][j];	
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
      } 
						
      
      
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
