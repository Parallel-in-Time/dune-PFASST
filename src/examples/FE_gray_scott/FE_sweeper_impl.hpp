
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
//#include <dune/functions/functionspacebases/interpolate.hh>

//#include <boost/math/constants/constants.hpp>

//using boost::math::constants::pi;
//using boost::math::constants::two_pi;
//using boost::math::constants::pi_sqr;

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
        /*config::options::add_option<size_t>("Heat FE", "num_dofs",
                                            "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat FE", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat FE", "nu",
                                               "thermal diffusivity");*/
      }


        template<class SweeperTrait, typename Enabled>
      Heat_FE<SweeperTrait, Enabled>::Heat_FE(std::shared_ptr<fe_manager> FinEl, size_t nlevel)
        :   IMEX<SweeperTrait, Enabled>()

      {

	
		    this->FinEl = FinEl;
	      basis = FinEl->get_basis(nlevel);
        grid  = FinEl->get_grid();

        assembleProblem(basis, this->A_dune, this->M_dune);

        const auto bs = basis->size();
        std::cout << "Finite Element basis consists of " <<  basis->size() << " elements " << std::endl;

        this->encap_factory()->set_size(bs);
        /*const unsigned int dim = 2;
        Dune::FieldVector<typename GridType::ctype,dim> L;
        L[0]=1; L[1]=1;
        typename Dune::array<int,dim> s;
        std::fill(s.begin(), s.end(), nelements);
        std::bitset<dim> periodic;//(true, true);
        periodic[0]=true; periodic[1]=true;

        grid        = std::make_shared<GridType>(L,s,periodic,0);


        grid->globalRefine(finer);
        this->finer = finer;

        typedef GridType::LeafGridView GridView;
        GridType::LeafGridView gridView = grid->leafGridView();

        std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
        std::cout << "***** Ordnung der Basis " << SweeperTrait::BASE_ORDER << std::endl;

        this->basis = std::make_shared<BasisFunction>(gridView);

        std::cout << "***** Basis erstellt mit " <<  basis->size() << " Elementen " << std::endl;

        this->encap_factory()->set_size(basis->size());*/
	
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::set_options()
      {


        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        int num_nodes = this->get_quadrature()->get_num_nodes();

        //assembleProblem(basis, A_dune, M_dune); //das muss bei mlsdc raus bei sdc rein!!


      }

        template<class SweeperTrait, typename Enabled>
        template<typename Basis>
        void
        Heat_FE<SweeperTrait, Enabled>::assemble(Basis &basis){
          //using namespace Dune;
          //typedef GridType::LeafGridView GridView;
          //using BasisFunction = Functions::PQkNodalBasis<GridView,1>; //SweeperTrait::BASE_ORDER>;
          //std::shared_ptr<BasisFunction> basis;
          //GridType::LeafGridView gridView       = grid->leafGridView();
          //this->basis       = std::make_shared<BasisFunction>(gridView);
          std::cout << "Im assemble mit basis" << std::endl;
          assembleProblem(basis, this->A_dune, this->M_dune);


        };

        template<class SweeperTrait, typename Enabled>
        void
        Heat_FE<SweeperTrait, Enabled>::assemble(){
          //using namespace Dune;
          //typedef GridType::LeafGridView GridView;
          //using BasisFunction = Functions::PQkNodalBasis<GridView,1>; //SweeperTrait::BASE_ORDER>;
          //std::shared_ptr<BasisFunction> basis;
          //GridType::LeafGridView gridView       = grid->leafGridView();
          //this->basis       = std::make_shared<BasisFunction>(gridView);
          std::cout << "Im assemble" << std::endl;
          basis;
          std::cout << "Im assemble nach basis aufruf" << std::endl;
          assembleProblem(basis, this->A_dune, this->M_dune);


        };



      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        const auto dim = SweeperTrait::DIM;
        spatial_t nu = this-> _nu; 
	
	
	auto exact_solution1 = [t,  nu, dim](const Dune::FieldVector<double,dim>&x){
             
	   
	   double eps = 1e-8;
	   
	  if( (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) < 0.0025)
	    return 0.5;
	    return 1.0;
            
        };

	auto exact_solution2 = [t, nu, dim](const Dune::FieldVector<double,dim>&x){
          double eps = 1e-8; 

	   if( (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) < 0.0025)
	   return 0.25;
	   return 0.0;
	};
	
	 auto N_x = [t](const Dune::FieldVector<double,dim>&x){
            return x;

        };

        Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
	Dune::BlockVector<Dune::FieldVector<double,1>> pom1, pom2;
	pom1.resize(basis->size());
	pom2.resize(basis->size());
        
	interpolate(*basis, x_node, N_x);

	interpolate(*basis, pom1, exact_solution1);
	interpolate(*basis, pom2, exact_solution2);
	
	
	
	for (int i = 0; i< basis->size(); ++i)
	{
	  result->data()[i][0] = pom1[i];
	  result->data()[i][1] = pom2[i];
	  
	  
	}

        //auto grid = (*(this->get_fine())).get_grid();
        typedef Dune::BlockVector<Dune::FieldVector<double, 2> > VectorType;
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ColumnType;
        typedef Dune::YaspGrid<2> GridType; //ruth_dim
        //typedef GridType::LevelGridView GridView;
        //GridType::LevelGridView gridView = grid->leafGridView();
        //Dune::VTKWriter<GridView> vtkWriter(gridView);



        /*string name = "initial";

        ColumnType sol_u, sol_v;
        sol_u.resize(result->data().size());
        sol_v.resize(result->data().size());
        for (int i =0; i< result->data().size(); ++i)
        {
          sol_u[i] = result->data()[i][0];
          sol_v[i] = result->data()[i][1];
        }


        vtkWriter.addVertexData(sol_u, "fe_solution_u");
        vtkWriter.addVertexData(sol_v, "fe_solution_v");


        vtkWriter.write("gray_scott" + name);*/




        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        this->_num_expl_f_evals = 0;
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
                     "  t["<<m<<"]=" << (t + dt * nodes[m])
                     << "      |abs residual| = " <<  this->_abs_res_norms[m]
                     << "      |rel residual| = " <<  this->_rel_res_norms[m]
//                      << "      |abs error| = " << encap::norm0(error[m])
//                      << "      |rel error| = " << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" << (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " << this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " << this->_rel_res_norms[num_nodes]
//                   << "      |abs error| = " << encap::norm0(error[num_nodes])
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

      
      typedef Dune::YaspGrid<2> GridType; //ruth_dim
      
      template<class SweeperTrait, typename Enabled>
      shared_ptr<GridType>
      Heat_FE<SweeperTrait, Enabled>::get_grid() const
      {
        return grid;
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
	    rel_error[m]=0;
        }

        return rel_error;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(), "evaluating EXPLICIT part at t=" << t);

        auto result = this->get_encap_factory().create();
	auto Mresult = this->get_encap_factory().create();
	
        result->zero();
	
	for(int i=0; i<u->data().size();++i)
	{
	  result->data()[i][0] = -u->data()[i][0]*u->data()[i][1]*u->data()[i][1] + this->_f*(1-u->data()[i][0]);	  	  
	  result->data()[i][1] =  u->data()[i][0]*u->data()[i][1]*u->data()[i][1] - (this->_f + this->_k )*u->data()[i][1];
  
	  
	}
	
	
	
	
	
        this->_num_expl_f_evals++;
	
	this->M_dune.mv(result->data(), Mresult->data());
	
	
	

	
        return Mresult;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {

        ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);


        auto result = this->get_encap_factory().create();

        double nu =this->_nu;

        this->A_dune.mmv(u->get_data(), result->data());



        /*for (size_t i = 0; i < result->data().size(); i++) {


          std::cout << "f->data im evaluate  " << rhs->data()[i] << std::endl;
        }*/
        //auto isDirichlet = [] (auto x) {return (x[0]<1e-8 or x[0]>0.9999 or x[1]<1e-8 or x[1]>0.9999);}; //ruth_dim


        /*Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dune);

        Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dune,1.0);

        Dune::CGSolver<VectorType> cg(linearOperator,
                                preconditioner,
                                1e-10, // desired residual reduction factor
                                500,    // maximum number of iterations
                                0);    // verbosity of the solver
        Dune::InverseOperatorResult statistics ;

        cg.apply(result->data(), rhs->data() , statistics );



        auto var = this->get_encap_factory().create();
        M_dune.mv(result->data(), var->data());
        auto neu = this->get_encap_factory().create();
        rhs2->scaled_add(-1.0 , var)  ;*/

        result->data() *= nu;

        /*for (size_t i = 0; i < result->data().size(); i++) {


          std::cout << "f->data2 im evaluate  " << result->data()[i] << std::endl;
        }*/

	/*for(int i=0; i<u->data().size();++i)
	{
	  std::cout << result->data()[i][0] << std::endl;	  	  
	  std::cout << result->data()[i][1] << std::endl;	  
  
	  
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
        Dune::BlockVector<Dune::FieldVector<double,SweeperTrait::NR_OF_COMP> > M_rhs_dune ;
        M_rhs_dune.resize(rhs->get_data().size());


	
        //M_dune.mv(rhs->data(), M_rhs_dune); //multipliziert rhs mit matrix_m_dune
		M_rhs_dune = rhs->get_data(); 
	

        Dune::BCRSMatrix<Dune::FieldMatrix<double,SweeperTrait::NR_OF_COMP,SweeperTrait::NR_OF_COMP> > M_dtA_dune = 	Dune::BCRSMatrix<Dune::FieldMatrix<double,2,2> >(this->A_dune);
        M_dtA_dune *= (dt /** this->_nu*/);
        M_dtA_dune += this->M_dune;
	
	
	

        Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dtA_dune);

        Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dtA_dune,1.0);

        Dune::BiCGSTABSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-15, // desired residual reduction factor
                              500,    // maximum number of iterations
                              0);    // verbosity of the solver
	
	
	typedef Dune::FieldVector<double, 2> Block;
        Dune::InverseOperatorResult statistics ;
	
	cg.apply(u->data(), M_rhs_dune , statistics ); //rhs ist nicht constant!!!!!!!!!

	
	
	std::cout << " " << std::endl; 
        ML_CVLOG(4, this->get_logger_id(),
                 "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);


	Dune::BlockVector<Dune::FieldVector<double,2> > M_u;
        M_u.resize(u->get_data().size());
	
	this->M_dune.mv(u->get_data(), M_u);
	
	/*std::cout << "********* mu *************************" << std::endl;
	for (int i=0; i< u->data().size(); i++){
	  std::cout <<  "u  " << u->data()[i] << std::endl;
          std::cout <<  "MU " << M_u[i] << std::endl;
	  std::cout <<  "M  " << this->M_dune[i][i] << std::endl;
        }
	std::cout << "**********************************" << std::endl;*/	
	
	//this->compute_residuals(false);
	//this->test_M(M_dune);
	//std::exit(0);
	
	
        for (size_t i = 0; i < u->data().size(); i++) {
	  Block b = (M_u[i] - rhs->data()[i]) ;
	  b/=dt;
          f->data()[i] = b;
          //std::cout << "f_solve " << f->data()[i] << std::endl;
        }
        
	/*std::cout << "#############################################################################" << std::endl;        
        for (size_t i = 0; i < u->get_data().size(); i++) {
		std::cout << u->data()[i][1] << std::endl; 
        } 
	std::cout << "#############################################################################" << std::endl;*/
        
        
        //evaluate_rhs_impl(0,u);
        //std::exit(0);
        this->_num_impl_solves++; 


	/*
	std::cout << "eval implicit "<<  std::endl;	
	evaluate_rhs_impl(t, u);
	std::exit(0);*/
	
	

      }
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
