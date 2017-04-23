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

//#include <leathers/push>
//#include <leathers/all>
#include <boost/math/constants/constants.hpp>
//#include <leathers/pop>
using boost::math::constants::pi;
using boost::math::constants::two_pi;
using boost::math::constants::pi_sqr;

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>

#include <iostream>
#include <c++/4.8/memory>

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
        /*config::options::add_option<size_t>("Heat FE", "num_dofs",
                                            "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat FE", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat FE", "nu",
                                               "thermal diffusivity");*/
      }


        template<class SweeperTrait, typename Enabled>
      Heat_FE<SweeperTrait, Enabled>::Heat_FE(std::shared_ptr<Dune::Functions::PQkNodalBasis<GridType::LeafGridView,SweeperTrait::BASE_ORDER>> basis, size_t nlevel)
        :   IMEX<SweeperTrait, Enabled>()

      {
      
	/*Dune::FieldVector<double,SweeperTrait::DIM> hR = {20};
	Dune::FieldVector<double,SweeperTrait::DIM> hL = {-20};
        array<int,SweeperTrait::DIM> n;

	std::fill(n.begin(), n.end(), nelements);

        this->grid = std::make_shared<GridType>(hL, hR, n);
        //grid.globalRefine(0);
	
        typedef GridType::LeafGridView GridView;
        GridType::LeafGridView gridView = grid->leafGridView();

        std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
        std::cout << "***** Ordnung der Basis " << SweeperTrait::BASE_ORDER << std::endl;
	
	
        this->basis = std::make_shared<BasisFunction>(gridView);

        std::cout << "***** Basis erstellt mit " <<  basis->size() << " Elementen " << std::endl;

        this->encap_factory()->set_size(basis->size());*/
	//this->FinEl = FinEl;
	//basis = FinEl->get_basis(nlevel);
	    
	    this->basis = basis;
	assembleProblem(basis, this->A_dune, this->M_dune);

        const auto bs = basis->size();
        std::cout << "Finite Element basis of level " << nlevel << " consists of " <<  basis->size() << " elements " << std::endl;

        this->encap_factory()->set_size(bs);

      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::set_options()
      {


        std::cout << " set_options " <<  std::endl;  
        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        int num_nodes = this->get_quadrature()->get_num_nodes();

        //assembleProblem(basis, A_dune, M_dune);
        std::cout << " set_options " <<  std::endl; 

      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        
                const auto dim = 1; //SweeperTrait::DIM;
	spatial_t n  = this-> _n;
    spatial_t l0 = this-> _nu;
	spatial_t l1 = l0/2. *(pow((1+n/2.), 1/2.) + pow((1+ n/2.), -1/2.) );
	spatial_t d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
	//std::cout << "nu = " << this->_nu << std::endl;
	//std::cout << "delta = " << this->_delta << std::endl;
        auto exact_solution = [l0, l1, n, d, t](const Dune::FieldVector<double,dim>&x){ 
	  return pow((1 + (pow(2, n/2.)-1 )* exp(-(n/2.)*d*(x+2*l1*t)) ), -2./n);
        };

        
        
        
        
        
        /*const auto dim = 1; //SweeperTrait::DIM;
        spatial_t nu = this-> _nu; 
	spatial_t delta = this->_delta;
	//std::cout << "nu = " << this->_nu << std::endl;
	//std::cout << "delta = " << this->_delta << std::endl;
        auto exact_solution = [t, nu, dim, delta](const Dune::FieldVector<double,dim>&x){
          double c = 2./delta;  
	  return 0.5*(1.0-std::tanh((x[0] - c*t*nu)/(delta)));
        };*/


	

        auto N_x = [t](const Dune::FieldVector<double,dim>&x){
            return x;

        };

        Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
        interpolate(*basis, x_node, N_x);

        interpolate(*basis, result->data(), exact_solution);

	/*for (int i=0; i< result->data().size(); i++){
	 std::cout << "result = " << result->data()[i] << std::endl;
	}*/
        return result;
      }
      
      
      
      
      //______________________________
      
      
      
      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::source(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        const auto dim = 1;
        spatial_t nu = this-> _nu; // nu je valjda desni rub
	//std::cout << "nu = " << this->_nu << std::endl;
        auto exact_solution_source = [t, nu, dim](const Dune::FieldVector<double,dim>&x){
            double solution=1.0;
            //for(int i=0; i<SweeperTrait::DIM; i++){solution *=x[i];}    //
            
	    std::cout << "PI = " << PI << std::endl;
	    std::exit(0);
	    return 2*t + PI*PI*(std::sin(PI*x[0])*x[0] + std::sin(PI*x[1])*x[1])  - 2*PI*(std::cos(PI*x[0]) + std::cos(PI*x[1])) + pow(std::sin(PI*x[0])*x[0] + std::sin(PI*x[1])*x[1] + t*t, 2) ;
        };


	
	

        auto N_x = [t](const Dune::FieldVector<double,dim>&x){
            return x;

        };

        Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
        interpolate(*basis, x_node, N_x);

        interpolate(*basis, result->data(), exact_solution_source);


        return result;
      }*/
      
      
      
      
      
      
      
      
      
      //_______________________________
      
      
      
      
      
      
      
      
      
      
      
      
      
      

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        //ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        //this->_num_expl_f_evals = 0;
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
//                      << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
//                      << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" <<  (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " <<  this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " <<  this->_rel_res_norms[num_nodes]
//                   << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
//                   << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes])
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

      //typedef Dune::YaspGrid<1,Dune::EquidistantOffsetCoordinates<double, 1> > GridType; //ruth_dim
     
     
      
      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<GridType>
      Heat_FE<SweeperTrait, Enabled>::get_grid() const
      {
        return grid;
      }*/


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

      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(),  "evaluating EXPLICIT part at t=" << t);

        auto result2 = this->get_encap_factory().create();
	auto Mresult = this->get_encap_factory().create();
        result2->zero();

	
	
        this->_num_expl_f_evals++;
	
        return result2;

      }*/

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
	result->data()*=_nu*_nu;
	this->A_dune.umv(u->get_data(), result->data());

	

        //result->data() *= nu;
	/*std::cout << "evaluate  " << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {

	  std::cout << "f " << result->data()[i] << std::endl;

         }*/
        
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
	Dune::BlockVector<Dune::FieldVector<double,1> > newton_rhs, newton_rhs2 ;
    newton_rhs.resize(rhs->get_data().size());
    newton_rhs2.resize(rhs->get_data().size());
    
	
	for (int i=0; i< 20 ;i++){
	  Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune); ///////M
	  evaluate_f(f, u, dt, rhs);
	  evaluate_df(df, u, dt);
	  df.mv(u->data(), newton_rhs);
	  newton_rhs -= f->data();
      newton_rhs2 = newton_rhs;

	  auto isLeftDirichlet = [] (auto x) {return (x[0] < -20.0 + 1e-8 ) ;};
	  auto isRightDirichlet = [] (auto x) {return (x[0] > 20.0 - 1e-8 ) ;};
 
	
	
	  std::vector<double> dirichletLeftNodes;
	  interpolate(*basis, dirichletLeftNodes, isLeftDirichlet);  

	  std::vector<double> dirichletRightNodes;
	  interpolate(*basis, dirichletRightNodes, isRightDirichlet);
	
	  /*for(int i=0; i<rhs->data().size(); ++i){
	
	    if(dirichletLeftNodes[i])
	    newton_rhs[i] = exact(t)->get_data()[i]; //1;//-dt;

	    if(dirichletRightNodes[i])
	    newton_rhs[i] = exact(t)->get_data()[i]; //0;


	  
	  }
	  for (size_t j=0; j<df.N(); j++){
	    if (dirichletRightNodes[j] || dirichletLeftNodes[j]){
	      auto cIt = df[j].begin();
	      auto cEndIt = df[j].end();
	      for(; cIt!=cEndIt; ++cIt){
		  *cIt = (j==cIt.index()) ? 1.0 : 0.0;
	      }
	    }
	  }*/
	  //std::cout << "5"<<std::endl;

  
	  Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df);
	  Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df,1.0);
	  Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-10, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver
	  Dune::InverseOperatorResult statistics ;
	  cg.apply(u->data(), newton_rhs , statistics ); //rhs ist nicht constant!!!!!!!!!
	  df.mv(u->data(), residuum->data());
      residuum->data() -= newton_rhs2;
      std::cout << "residuums norm " << residuum->norm0() << std::endl;
      if (residuum->norm0()< _abs_newton_tol){break;}
	  
        for (size_t i = 0; i < u->get_data().size(); i++) {

	  //std::cout << "u " << u->data()[i] << std::endl;

        }
	}

	
	
	//std::exit(0);
	
        /*for (size_t i = 0; i < u->get_data().size(); i++) {

	  std::cout << "u " << u->data()[i] << std::endl;

        }*/
	
	
	
	
	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);

//std::cout << "impl solve "  << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
	  //std::cout << "f " << f->data()[i] << std::endl;
        }
        //evaluate_rhs_impl(0, u);
	//std::exit(0);
        this->_num_impl_solves++;
        //if (this->_num_impl_solves==5) std::exit(0);


      }
      
            template<class SweeperTrait, typename Enabled>
      void
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
	  fneu[i]= pow(u->get_data()[i], _n+1) - u->get_data()[i];	
	}
	f->data() *= (_nu*_nu);
	this->M_dune.mv(fneu, f->data());
	this->A_dune.mmv(u->get_data(),f->data());
	f->data() *= dt;
	this->M_dune.umv(u->get_data(),f->data());
	f->data() -=rhs->get_data();
          
	
	/*f->zero();
	
	Dune::BlockVector<Dune::FieldVector<double,1> > fneu;
        fneu.resize(u->get_data().size());
	for (int i=0; i<f->data().size(); ++i)
	{fneu[i]= -8*this->_nu*this->_nu *u->data()[i]*u->data()[i]*(1.00-u->data()[i])/(this->_delta*this->_delta);

	}
	
	this->M_dune.mv(fneu, f->data());
	
	this->A_dune.mmv(u->get_data(),f->data());
	f->data() *= dt;
	this->M_dune.umv(u->get_data(),f->data());
	f->data() -=rhs->get_data();*/

      }
						
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						){
          
          
          
          
          
          for (int i=0; i<df.N(); ++i)
            {
            for (int j=0; j<df.M(); ++j)
                {
                    if (df.exists(i,j)) 
                        df[i][j]= (_nu*_nu)*(_n+1) *this->M_dune[i][j] * pow(u->get_data()[j], _n);	
                }
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
          
          
          
	/*for (int i=0; i<df.N(); ++i)
	{
	  for (int j=0; j<df.M(); ++j)
	    {
	  if (df.exists(i,j)) {
	    df[i][j]=  -(_nu*_nu)*this->M_dune[i][j]  *(16* ((double) u->get_data()[j]) -24* ((double) u->get_data()[j])*((double) u->get_data()[j]))/(this->_delta*this->_delta);	
	  }
	    //std::cout << df[i][j]<<std::endl;
	  }
	}
	df-=this->A_dune;
	
	df*=dt;
	df+=this->M_dune;*/

      }  
      
      
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
