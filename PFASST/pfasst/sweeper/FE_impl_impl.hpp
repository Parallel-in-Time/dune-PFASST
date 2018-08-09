//#include "pfasst/sweeper/FE_imex.hpp"

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <memory>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Core>

#include <ios>



#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/densematrix.hh>
#include<dune/common/parallel/mpihelper.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/io.hh>
#include<dune/istl/matrixmarket.hh>
#include<dune/istl/matrixredistribute.hh>
#include <dune/istl/schwarz.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/typetree/utility.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#if USE_DG
#  include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#else
#  include <dune/functions/functionspacebases/pqknodalbasis.hh>
#endif
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/istlbackend.hh>
#include <dune/fufem/formatstring.hh>
#include <dune/fufem/assemblers/transferoperatorassembler.hh>

#include<dune/istl/paamg/pinfo.hh>
#include<dune/istl/paamg/graph.hh>
/*#if USE_DG
#  include <dune/parmg/test/dglaplacematrix.hh>
#else
#  include <dune/parmg/test/laplacematrix.hh>
#endif*/
/*#include <dune/parmg/iterationstep/lambdastep.hh>
#include <dune/parmg/iterationstep/multigrid.hh>
#include <dune/parmg/iterationstep/multigridstep.hh>
#include <dune/parmg/norms/normadapter.hh>
#include <dune/parmg/parallel/communicationp1.hh>
#include <dune/parmg/parallel/communicationdg.hh>
#include <dune/parmg/parallel/datahandle.hh>
#include <dune/parmg/parallel/dofmap.hh>
#include <dune/parmg/parallel/globaldofindex.hh>
#include <dune/parmg/parallel/istlcommunication.hh>
#include <dune/parmg/parallel/matrixalgebra.hh>
#include <dune/parmg/parallel/vectoralgebra.hh>
#include <dune/parmg/parallel/redistributematrix.hh>
#include <dune/parmg/parallel/redistributevector.hh>
#include <dune/parmg/parallel/parallelenergyfunctional.hh>
#include <dune/parmg/parallel/parallelenergynorm.hh>
#include <dune/parmg/parallel/restrictmatrix.hh>
#include <dune/parmg/solvers/coarsesuperlusolver.hh>
#include <dune/parmg/solvers/linesearch.hh>
#include <dune/parmg/solvers/directionsearch.hh>
#include <dune/parmg/iterationstep/multigridsetup.hh>
#include <dune/parmg/iterationstep/parallelprojectedgs.hh>

#include <dune/solvers/iterationsteps/blockgssteps.hh>
#include <dune/solvers/norms/energynorm.hh>
#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/transferoperators/compressedmultigridtransfer.hh>

using namespace Dune::ParMG;*/

namespace pfasst
{
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  IMEX<SweeperTrait, BaseFunction, Enabled>::IMEX()
    :   Sweeper<SweeperTrait, BaseFunction, Enabled>()
      , _q_integrals(0)
      , _impl_rhs(0)
      , _impl_rhs_restrict(0)
      , _num_impl_f_evals(0)
      , _num_impl_solves(0)
  {}

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::initialize()
  {
    pfasst::Sweeper<SweeperTrait, BaseFunction, Enabled>::initialize();
    //std::cout << "im initialize "<< std::endl;
    const auto num_nodes = this->get_quadrature()->get_num_nodes();
    //std::cout << this->get_states().size() << std::endl;
    //std::cout << num_nodes+1 << std::endl;
    assert(this->get_states().size() == num_nodes + 1);

    this->_q_integrals.resize(num_nodes + 1);
    std::generate(this->_q_integrals.begin(), this->_q_integrals.end(),
             std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

    this->_impl_rhs.resize(num_nodes + 1);
    std::generate(this->_impl_rhs.begin(), this->_impl_rhs.end(),
             std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));
    
    this->_impl_rhs_restrict.resize(num_nodes + 1);
    std::generate(this->_impl_rhs_restrict.begin(), this->_impl_rhs_restrict.end(),
             std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));
    

    this->compute_delta_matrices();
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::setup()
  {
    pfasst::Sweeper<SweeperTrait, BaseFunction, Enabled>::setup();

    ML_CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id(),
      "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.");
    //std::cout << "wird benutzt " << std::endl; std::exit(0);
    this->initialize();
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::pre_predict()
  {
    Sweeper<SweeperTrait, BaseFunction, Enabled>::pre_predict();

  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::predict()
  {
    

    Sweeper<SweeperTrait, BaseFunction, Enabled>::predict();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    ML_CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id(),
      "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.");

    const typename traits::time_t t = this->get_status()->get_time();
    const typename traits::time_t dt = this->get_status()->get_dt();

    assert(this->get_quadrature() != nullptr);
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), typename traits::time_t(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    this->_impl_rhs.front() = this->evaluate_rhs_impl(t, this->get_states().front());


    ML_CLOG(INFO, this->get_logger_id(),  "Predicting from t=" << t << " over " << num_nodes << " nodes"
                          << " to t=" << (t + dt) << " (dt=" << dt << ")");
    typename traits::time_t tm = t;
    	//int my_rank, num_pro;
        //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        //MPI_Comm_size(MPI_COMM_WORLD, &num_pro );MPI_Barrier(MPI_COMM_WORLD); 
    for (size_t m = 0; m < num_nodes; ++m) {
      this->states()[m + 1]->data() = this->states()[m]->data();
    //if(my_rank==0) for(auto r =this->states()[m]->data().begin(); r !=this->states()[m]->data().end(); ++r){std::cout << "ds states " << *r <<std::endl;} //std::exit(0);
      //std::cout << "vor zuweisung " << std::endl; //std::exit(0);	
      this->_impl_rhs[m + 1] = this->evaluate_rhs_impl(tm, this->get_states()[m + 1]);
      //std::cout << "nach zuweisung " << std::endl; std::exit(0);	
      tm += dt *  (nodes[m+1] - nodes[m]);
      ML_CVLOG(1, this->get_logger_id(), "");

    }       //for(auto r =this->_impl_rhs[2]->data().begin(); r !=this->_impl_rhs[2]->data().end(); ++r){std::cout << "predict rhs_impl " << *r <<std::endl;} std::cout << "num_nodes " << num_nodes <<std::endl; //std::exit(0);	

  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::post_predict()
  {
    //for(auto r =this->_impl_rhs[2]->data().begin(); r !=this->_impl_rhs[2]->data().end(); ++r){std::cout << "postpredict rhs_impl " << *r <<std::endl;} //std::exit(0);	
    Sweeper<SweeperTrait, BaseFunction, Enabled>::post_predict();
    //for(auto r =this->_impl_rhs[2]->data().begin(); r !=this->_impl_rhs[2]->data().end(); ++r){std::cout << "postpredict2 rhs_impl " << *r <<std::endl;} //std::exit(0);	
  }

  
  
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::pre_sweep()
  {
    
    //for(auto r =this->_impl_rhs[2]->data().begin(); r !=this->_impl_rhs[2]->data().end(); ++r){std::cout << "presweep rhs_impl " << *r <<std::endl;} std::exit(0);	
    Sweeper<SweeperTrait, BaseFunction, Enabled>::pre_sweep();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    ML_CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id(),
      "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.");



	
    const typename traits::time_t dt = this->get_status()->get_dt();
    const auto q_mat = this->get_quadrature()->get_q_mat();
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), typename traits::time_t(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    ML_CVLOG(4, this->get_logger_id(), "computing integrals");
    ML_CVLOG(6, this->get_logger_id(), "  q_int     = dt * Q * f_ex");
    //this->_q_integrals = encap::mat_mul_vec(dt, q_mat, this->_expl_rhs);
    ML_CVLOG(6, this->get_logger_id(), "           += dt * Q * f_im");

    	int my_rank, num_pro;
        //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        //MPI_Comm_size(MPI_COMM_WORLD, &num_pro );MPI_Barrier(MPI_COMM_WORLD); 
    //if(my_rank==0) for(auto r =this->states()[2]->data().begin(); r !=this->states()[2]->data().end(); ++r){std::cout << "ds states " << *r <<std::endl;}// std::exit(0);
    //for(auto r =this->_impl_rhs[2]->data().begin(); r !=this->_impl_rhs[2]->data().end(); ++r){std::cout << "pre_sweep rhs_impl " << *r <<std::endl;} //std::exit(0);

    encap::mat_apply(this->_q_integrals, dt, q_mat, this->_impl_rhs, true); //false

    ML_CVLOG(4, this->get_logger_id(), "  subtracting function evaluations of previous iteration and adding FAS correction");


    for (size_t m = 0; m < num_nodes; ++m) {
      for (size_t n = 0; n < m + 1; ++n) {
	//for(auto r =this->_impl_rhs[n + 1]->data().begin(); r !=this->_impl_rhs[n + 1]->data().end(); ++r){std::cout << "M result " << *r <<std::endl;} std::exit(0);
        this->_q_integrals[m + 1]->scaled_add(-dt * this->_q_delta_impl(m + 1, n + 1), this->_impl_rhs[n + 1]);
      }

//       ML_CVLOG(6, this->get_logger_id(), LOG_FLOAT << "  q_int["<<(m+1)<<"] += tau["<<(m+1)<<"]                  = "
//                                        << to_string(this->get_tau()[m + 1]));
      this->_q_integrals[m + 1]->scaled_add(1.0, this->get_tau()[m + 1]);
    }

//     for (size_t m = 0; m < num_nodes + 1; ++m) {
//       ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "  q_int["<<m<<"] = " << to_string(this->_q_integrals[m]));
//     }
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::sweep()
  {
    
    
    Sweeper<SweeperTrait, BaseFunction, Enabled>::sweep();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    ML_CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id(),
      "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.");

    const typename traits::time_t t = this->get_status()->get_time();
    const typename traits::time_t dt = this->get_status()->get_dt();
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), typename traits::time_t(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    ML_CLOG(INFO, this->get_logger_id(), "Sweeping from t=" << t << " over " << num_nodes
                          << " nodes to t=" << (t + dt) << " (dt=" << dt << ")");

    //this->_expl_rhs.front() = this->evaluate_rhs_expl(t, this->get_states().front());

    typename traits::time_t tm = t;
    // note: m=0 is initial value and not a quadrature node
    for (size_t m = 0; m < num_nodes; ++m) {
      ML_CVLOG(4, this->get_logger_id(), "propagating from t["<<m<<"]=" << (t + (dt * nodes[m]))
                                                   << " to t["<<(m+1)<<"]=" << (t + (dt * nodes[m+1])));
//       ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "  u["<<m<<"] = " << to_string(this->get_states()[m]));

      // compute right hand side for implicit solve (i.e. the explicit part of the propagation)
      //std::cout << "vorm create rhs" << std::endl;	
      shared_ptr<typename traits::encap_t> rhs = this->get_encap_factory().create();
      std::cout << "nach create rhs " << std::endl;
      // rhs = u_0
      
      if (is_coarse){
        rhs->data() =  this->_M_initial->get_data(); //   this->get_states().front()->get_data();	
	
	//for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << "sweeper " << *r <<std::endl;}
	//MPI_Barrier(MPI_COMM_WORLD);
	//std::cout << "rhs debugging " << std::endl;
	//std::exit(0);	

          
      }else{
        //M_dune.mv(this->get_states().front()->get_data(), rhs->data()); //hier ersetze ich jetzt duch allgemeines apply_Mass ruth
	//for(auto r =this->get_states().front()->data().begin(); r != this->get_states().front()->data().end(); ++r){std::cout << "rhs vorm mult " << *r <<std::endl;}
	//std::cout << "rhs debugging " << std::endl;
	//M_dune->jacobian_apply((this->get_states().front()->data()), (rhs->data()));
	//std::cout << "rhs glueck " << std::endl;
	this->get_states().front()->apply_Mass(M_dune, rhs); 	//std::cout << "bevor ich den quatsch hier lese" << std::endl;
	//for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << "rhs nachm mult " << *r <<std::endl;}
	//std::cout << "rhs debugging " << std::endl;
	//std::exit(0);


	//this->get_states().front()->apply_Mass(M_dune);
        

        
      }
      //std::cout << "is coarse " << is_coarse <<  std::endl; std::exit(0);
      
      //ML_CVLOG(6, this->get_logger_id(), "  rhs = u[0]                    = " << to_string(rhs));

      // rhs += dt * \sum_{i=0}^m (QI_{m+1,i} fI(u_i^{k+1}) + QE_{m+1,i-1} fE(u_{i-1}^{k+1}) ) + QE_{m+1,m} fE(u_{m}^{k+1})
      for (size_t n = 0; n <= m; ++n) {
        rhs->scaled_add(dt * this->_q_delta_impl(m + 1, n), this->_impl_rhs[n]);
	//std::cout << dt * this->_q_delta_impl(m + 1, n) << std::endl;
      } //       MPI_Barrier(MPI_COMM_WORLD);// std::exit(0);

      //for(auto r =this->_q_integrals[m + 1]->data().begin(); r !=this->_q_integrals[m + 1]->data().end(); ++r){std::cout << "q_integrals " << *r <<std::endl;}
      //for(auto r =this->_q_integrals[m + 1]->data().begin(); r != this->_q_integrals[m + 1]->data().end(); ++r){std::cout << "q_integrals " << *r <<std::endl;}
      rhs->scaled_add(1.0, this->_q_integrals[m + 1]);
    	int my_rank, num_pro;
        //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        //MPI_Comm_size(MPI_COMM_WORLD, &num_pro );MPI_Barrier(MPI_COMM_WORLD); 
	//if (my_rank==0) for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << "rhs nachm mult " << *r <<std::endl;}     MPI_Barrier(MPI_COMM_WORLD); 
	//if (my_rank==1) for(auto r =rhs->data().begin(); r !=rhs->data().end(); ++r){std::cout << "rhs nachm mult " << *r <<std::endl;}     MPI_Barrier(MPI_COMM_WORLD); std::exit(0);
      //std::exit(0);
      // solve the implicit part
      ML_CVLOG(4, this->get_logger_id(), "  solve(u["<<(m+1)<<"] - dt * QI_{"<<(m+1)<<","<<(m+1)<<"} * f_im["<<(m+1)<<"] = rhs)");
      std::cout << "vor solve " << std::endl;
      this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, dt * this->_q_delta_impl(m+1, m+1), rhs);
      std::cout << "nach solve" << std::endl;
    //std::cout << "ende aufruuf impl solve " << std::endl;
      // reevaluate the explicit part with the new solution value
      tm += dt * this->_q_delta_impl(m+1, m+1);
      
      ML_CVLOG(4, this->get_logger_id(), "");
    //std::cout << "schleife " << num_nodes << " " << m << std::endl;

    }
    //std::cout << "ende sweep " << std::endl;
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::post_sweep()
  {
    Sweeper<SweeperTrait, BaseFunction, Enabled>::post_sweep();
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::post_step()
  {
    Sweeper<SweeperTrait, BaseFunction, Enabled>::post_step();
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::advance(const size_t& num_steps)
  {
    UNUSED(num_steps);

    assert(this->get_end_state() != nullptr);
    ML_CVLOG(4, this->get_logger_id(), "advancing");


    this->initial_state()->data() = this->get_end_state()->get_data();
    assert(this->get_quadrature() != nullptr);

    if (this->get_quadrature()->left_is_node() && this->get_quadrature()->right_is_node()) {
      //assert(this->_expl_rhs.front() != nullptr && this->_expl_rhs.back() != nullptr);
      assert(this->_impl_rhs.front() != nullptr && this->_impl_rhs.back() != nullptr);

      //this->_expl_rhs.front()->data() = this->_expl_rhs.back()->get_data();
      this->_impl_rhs.front()->data() = this->_impl_rhs.back()->get_data();
    } else {
      // TODO: this might not be necessary as it is probably dealt with in pre_predict and pre_sweep
//       throw NotImplementedYet("advancing IMEX for nodes not containing left and right time interval borders");
    }
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::advance()
  {
    this->advance(1);
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::reevaluate(const bool initial_only)
  {
    assert(this->get_status() != nullptr);
    assert(this->get_quadrature() != nullptr);

    const typename traits::time_t t0 = this->get_status()->get_time();

    if (initial_only) {
      assert( this->_impl_rhs.front() != nullptr);
      //this->_expl_rhs.front() = this->evaluate_rhs_expl(t0, this->get_initial_state());
	//std::cout << this->get_initial_state()[0] << std::endl;
	//std::exit(0);
        this->_impl_rhs.front() = this->evaluate_rhs_impl(t0, this->get_initial_state());

    } else {
      const typename traits::time_t dt = this->get_status()->get_dt();
      auto nodes = this->get_quadrature()->get_nodes();
      nodes.insert(nodes.begin(), 0.0);
	//std::exit(0);
      for (size_t m = 0; m < this->get_quadrature()->get_num_nodes() + 1; ++m) {
        const typename traits::time_t t = t0 + dt * nodes[m];
        assert( this->_impl_rhs[m] != nullptr);

        //this->_expl_rhs[m] = this->evaluate_rhs_expl(t, this->get_states()[m]);
        this->_impl_rhs[m] = this->evaluate_rhs_impl(t, this->get_states()[m]);
      }
    }
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::reevaluate()
  {
    this->reevaluate(false);
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>
  IMEX<SweeperTrait, BaseFunction, Enabled>::integrate(const typename SweeperTrait::time_t& dt)
  {
    auto const q_mat = this->get_quadrature()->get_q_mat();

    //auto result = encap::mat_mul_vec(dt, q_mat, this->_expl_rhs);
    auto result = encap::mat_mul_vec(dt, q_mat, this->_impl_rhs);
    //encap::mat_apply(result, dt, q_mat, this->_impl_rhs, false);

    
    //std::cout << "groesse von rhs_impl " << this->_impl_rhs[0]->get_data().size() << std::endl;
    return result;
  }
  
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>
  IMEX<SweeperTrait, BaseFunction, Enabled>::integrate_new(const typename SweeperTrait::time_t& dt)
  {
    auto const q_mat = this->get_quadrature()->get_q_mat();

    //vector<shared_ptr<typename SweeperTrait::encap_t>>&
    
    auto result = encap::mat_mul_vec(dt, q_mat, this->_impl_rhs_restrict);
    
    
    //std::cout << "groesse von rhs_impl " << this->_impl_rhs[0]->get_data().size() << std::endl;
    return result;
  }
  
  


  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::integrate_end_state(const typename SweeperTrait::time_t& dt)
  {
    try {
      Sweeper<SweeperTrait, BaseFunction, Enabled>::integrate_end_state(dt);
    } catch (std::runtime_error err) {
      assert(this->get_quadrature() != nullptr);
      assert(this->get_initial_state() != nullptr);

      this->end_state()->data() = this->get_initial_state()->get_data(); //achtung
      //M_dune.mv(this->get_initial_state()->get_data(), this->end_state()->data());
      //this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_expl_rhs)[0]);
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_impl_rhs)[0]);
//       ML_CVLOG(1, this->get_logger_id(), "end state: " << to_string(this->get_end_state()));
    }
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction,  Enabled>::compute_residuals(const bool& only_last)
  {
    ML_CVLOG(4, this->get_logger_id(), "computing residuals");

    assert(this->get_status() != nullptr);
    assert(this->get_quadrature() != nullptr);
    assert(this->get_initial_state() != nullptr);
    //int rank, num_pro;
    //MPI_Comm_rank(this->comm, &rank );
    const typename traits::time_t dt = this->get_status()->get_dt();
    const size_t num_nodes = this->get_quadrature()->get_num_nodes() + 1;
    //if (rank==0) std::cout << "***************************************************************  compute residuuals " << std::endl;
    if (only_last) {
      //if (rank==0) std::cout << "***************************************************************  compute residuuals only last " << std::endl;
      const size_t cols = this->get_quadrature()->get_q_mat().cols();
      const size_t rows = this->get_quadrature()->get_q_mat().rows();
      
      
      
      if (is_coarse){
        //if (rank==0) std::cout << "***************************************************************  compute residuuals only last coarse" << std::endl;
        this->residuals().back()->data() =  this->_M_initial->get_data(); //   this->get_states().front()->get_data();
      }else{
        //if (rank==0) std::cout << "***************************************************************  compute residuuals only last fine" << std::endl;
        //M_dune.mv(this->get_initial_state()->get_data(), this->residuals().back()->data());
	this->get_initial_state()->apply_Mass(M_dune, this->residuals().back());
      }
      
      //M_dune.mv(this->get_initial_state()->get_data(), this->residuals().back()->data());
      //this->residuals().back()->data() = this->get_initial_state()->get_data();
      
      shared_ptr<typename traits::encap_t> uM = this->get_encap_factory().create();	

      //M_dune.mv(this->get_states().back()->get_data(), uM->data());
      this->get_states().back()->apply_Mass(M_dune,uM);	
      //this->get_states().back()->apply_Mass();
      this->get_states().back()->apply_Mass(M_dune, uM);	
      this->residuals().back()->scaled_add(-1.0, uM);	
      //this->residuals()[m]->scaled_add(-1.0,uM);

      
      //for(auto r =this->_impl_rhs[2]->data().begin(); r !=this->_impl_rhs[2]->data().end(); ++r){std::cout << "residuals1 rhs_impl " << *r <<std::endl;} //std::exit(0);
      this->residuals().back()->scaled_add(1.0, this->get_tau().back());
      for (size_t n = 0; n < cols; ++n) {
        //this->residuals().back()->scaled_add(dt * this->get_quadrature()->get_q_mat()(rows - 1, n), this->_expl_rhs[n]);
        this->residuals().back()->scaled_add(dt * this->get_quadrature()->get_q_mat()(rows - 1, n), this->_impl_rhs[n]);
      }
      //for(auto r =this->_impl_rhs[2]->data().begin(); r !=this->_impl_rhs[2]->data().end(); ++r){std::cout << "residuals2 rhs_impl " << *r <<std::endl;} //std::exit(0);	
    } else {
      //if (rank==0) std::cout << "***************************************************************  compute residuuals NOT only last " << std::endl;

      for (size_t m = 0; m < num_nodes; ++m) {
        assert(this->get_states()[m] != nullptr);
        assert(this->residuals()[m] != nullptr);
        
        /*if(!is_coarse){
        for(int i=0; i< this->get_residuals()[m]->get_data().size(); i++){
           if(rank==0) std::cout<< rank << " rein ini " << m << " " << i << " " << this->get_initial_state()->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);
        for(int i=0; i< this->get_residuals()[m]->get_data().size(); i++){
           if(rank==1) std::cout<< rank << " rein ini " << m << " " << i << " " << this->get_initial_state()->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);}*/

    
    	if (is_coarse){
    		//std::cout << "------------------------------------------------------is coarse " << std::endl; 
        	this->residuals()[m]->data() =  this->_M_initial->get_data(); 
    	}else{
    		//std::cout << "------------------------------------------------------is fine " << std::endl; 
		this->get_initial_state()->apply_Mass(M_dune,  this->residuals()[m]);
    	}
    	/*MPI_Barrier(MPI_COMM_WORLD);
    	
    	for (int i=0; i<M_dune.N(); ++i)
            {
            for (int j=0; j<M_dune.M(); ++j)
                {
                    if(rank==0) if (M_dune.exists(i,j)) {std::cout << M_dune[i][j] << " ";}
                    //else {std::cout << std::endl;} 
                } if(rank==0) std::cout << std::endl;
            }       
    	
    	
    	MPI_Barrier(MPI_COMM_WORLD);
    	for (int i=0; i<M_dune.N(); ++i)
            {
            for (int j=0; j<M_dune.M(); ++j)
                {
                    if(rank==1) if (M_dune.exists(i,j)) {std::cout << M_dune[i][j] << " ";}
                    //else {std::cout << std::endl;} 
                }if(rank==1) std::cout << std::endl;
            }       
    	
    	
    	MPI_Barrier(MPI_COMM_WORLD);*/
    

        //if (!is_coarse) std::exit(0);    
    
	shared_ptr<typename traits::encap_t> uM = this->get_encap_factory().create();
	this->get_states()[m]->apply_Mass(M_dune, uM);
	this->residuals()[m]->scaled_add(-1.0,uM);
        assert(this->get_tau()[m] != nullptr);
        
        this->residuals()[m]->scaled_add(1.0, this->get_tau()[m]);
        
      }

      ML_CVLOG(5, this->get_logger_id(), "vor  res += dt * Q * F_im");

      /*for (size_t m = 0; m < num_nodes; ++m) {
        double d1 = this->get_residuals()[m]->norm0();
        double d2 =this->_impl_rhs[m]->norm0();
        if(rank==0){ ML_CVLOG(0, this->get_logger_id(), "  vor mat_apply  |res["<<m<<"]| = " << d1 << " "<< d2);}
        for(int i=0; i< this->get_residuals()[m]->get_data().size(); i++){
           if(rank==0) std::cout<< "get residuals " << m << " " << i << " " << this->get_residuals()[m]->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);
        for(int i=0; i< this->get_residuals()[m]->get_data().size(); i++){
           if(rank==1) std::cout<< "get residuals " << m << " " << i << " " << this->get_residuals()[m]->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);

      }MPI_Barrier(MPI_COMM_WORLD);*/

      encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_impl_rhs, false); //DAS IST FALSCH RUTH!!! 
      	
      	/*for (size_t m = 0; m < num_nodes; ++m) {
        for(int i=0; i< this->get_residuals()[m]->get_data().size(); i++){
           if(rank==0) std::cout<< rank << " residuals " << m << " " << i << " " << this->get_residuals()[m]->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);
        for(int i=0; i< this->get_residuals()[m]->get_data().size(); i++){
           if(rank==1) std::cout<< rank << " residuals " << m << " " << i << " " << this->residuals()[m]->get_data()[i] <<std::endl;         	
        }MPI_Barrier(MPI_COMM_WORLD);
        }*/
      	
      /*for (size_t m = 0; m < num_nodes; ++m) {
	if (rank!=num_pro-1) this->get_residuals()[m]->data()[this->get_residuals()[m]->data().size()-1] =0;
	if (rank!=0) this->get_residuals()[m]->data()[0] =0; 
        //double d1 = this->get_residuals()[m]->norm0();
        double d1_ = this->get_residuals()[m]->norm0(true, this->comm);
        //double d2 =this->_impl_rhs[m]->norm0();
        //double d2_ =this->_impl_rhs[m]->norm0(true, this->comm);
        //if(rank==0){ ML_CVLOG(0, this->get_logger_id(), " mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm nach mat_apply  |res["<<m<<"]| = " << d1_ << " ");}
        //typedef Dune::BlockVector<Dune::FieldVector<double,NR_COMP> > VectorType;
        //VectorType a = this->get_residuals()[m]->data();
        //collect(a);
        //if(rank==0){ ML_CVLOG(0, this->get_logger_id(), "  nach collect  |res["<<m<<"]| = " << d1 << " "<< d2);}
      }//MPI_Barrier(MPI_COMM_WORLD);*/      
    }
    //std::exit(0);
  }

  /*template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::test_m(MatrixType M)
  {
    std::cout << "************* im fe_imex pointer *********************" << std::endl;	
	for (int i=0; i< 2; i++){
	  //std::cout <<  "A  " << this->A_dune[i][i] << std::endl;
	  std::cout <<  "M  " << M[i][i] << std::endl;
        }
	std::cout << "**********************************" << std::endl;	
  }*/
  
  
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::compute_residuals()
  {
    this->compute_residuals(false);
  }

  /**
   * @throws std::runtime_error if not overwritten in specialized implementation
   */
  /*template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_t>
  IMEX<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u)
  {
    UNUSED(t); UNUSED(u);
    throw std::runtime_error("evaluation of explicit part of right-hand-side");
  }*/

  /**
   * @throws std::runtime_error if not overwritten in specialized implementation
   */

  

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_t>
  IMEX<SweeperTrait, BaseFunction, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u)
  {
    UNUSED(t); UNUSED(u);
    throw std::runtime_error("evaluation of implicit part of right-hand-side");
  }


  /*template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_t>
  IMEX<SweeperTrait, BaseFunction, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u)
  {
    UNUSED(t); UNUSED(u);
    throw std::runtime_error("evaluation of implicit part of right-hand-side");
  }*/

  /**
   * @throws std::runtime_error if not overwritten in specialized implementation
   */
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                              shared_ptr<typename SweeperTrait::encap_t> u,
                                              const typename SweeperTrait::time_t& t,
                                              const typename SweeperTrait::time_t& dt,
                                              const shared_ptr<typename SweeperTrait::encap_t> rhs)
  {
    UNUSED(f); UNUSED(u); UNUSED(t); UNUSED(dt); UNUSED(dt); UNUSED(rhs);
    throw std::runtime_error("spatial solver");
  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::compute_delta_matrices()
  {
    assert(this->get_quadrature() != nullptr);
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();
    auto nodes = this->get_quadrature()->get_nodes();
    if (this->get_quadrature()->left_is_node()) {
      ML_CLOG(ERROR, this->get_logger_id(),
              "Don't know how to compute delta matrices for quadrature containing left time point.");
      throw std::runtime_error("IMEX with quadrature containing left time point");
    } else {
      nodes.insert(nodes.begin(), typename traits::time_t(0.0));
    }

    ML_CVLOG(4, this->get_logger_id(), "computing Q_delta matrices for IMEX scheme");

    this->_q_delta_expl = Matrix<typename traits::time_t>::Zero(num_nodes + 1, num_nodes + 1);
    this->_q_delta_impl = Matrix<typename traits::time_t>::Zero(num_nodes + 1, num_nodes + 1);

    bool lr=true;
    if(lr){
        auto const q_mat = this->get_quadrature()->get_q_mat();



           Matrix<typename traits::time_t> q_mat_part(num_nodes , num_nodes );
           for (size_t m = 0; m < num_nodes ; ++m) {
             for (size_t n = 0; n < num_nodes ; ++n) {
               q_mat_part(m,n)=q_mat(m+1,n+1);
             }
           }



           Matrix<typename traits::time_t> q_mat_part_transpose=q_mat_part.transpose();



           Eigen::HouseholderQR<Matrix<typename traits::time_t>> qr(q_mat_part_transpose);
           Matrix<typename traits::time_t> qr_mat  = qr.householderQ();
           Matrix<typename traits::time_t> qr_matrix  = qr.matrixQR();

           Eigen::PartialPivLU<Matrix<typename traits::time_t>> lu(q_mat_part_transpose);
           Matrix<typename traits::time_t> u = lu.matrixLU();
           Matrix<typename traits::time_t> l = lu.matrixLU();


           for (size_t m = 0; m < num_nodes ; ++m) {
             for (size_t n = m; n < num_nodes; ++n) {
               l(m,n) = 0;
             }
           }
           for (size_t m = 0; m < num_nodes ; ++m) {
             l(m,m) = 1;
           }

           for (size_t m = 0; m < num_nodes ; ++m) {
             for (size_t n = 0; n < m; ++n) {
               qr_matrix(m,n) = 0;
             }
           }


           Matrix<typename traits::time_t> ut = u.transpose();



           for (size_t m = 1; m < num_nodes + 1; ++m) {
             for (size_t n = m; n < num_nodes + 1; ++n) {
               this->_q_delta_impl(n, m ) = ut(n-1, m-1);
             }
           }




           /*for (size_t m = 1; m < num_nodes + 1; ++m) {
             for (size_t n = m; n < num_nodes + 1; ++n) {
               this->_q_delta_expl(n, m - 1) = nodes[m] - nodes[m - 1];
               //this->_q_delta_impl(n, m) = nodes[m] - nodes[m - 1];
             }
           }*/



    }else{



        for (size_t m = 1; m < num_nodes + 1; ++m) {
            for (size_t n = m; n < num_nodes + 1; ++n) {
                //this->_q_delta_expl(n, m - 1) = nodes[m] - nodes[m - 1];
                this->_q_delta_impl(n, m) = nodes[m] - nodes[m - 1];
            }
        }


    }

    ML_CVLOG(5, this->get_logger_id(), "QE:");
    for (int row = 0; row < this->_q_delta_expl.rows(); ++row) {
      ML_CVLOG(5, this->get_logger_id(),
               "  " << this->_q_delta_expl.block(row, 0, 1, this->_q_delta_expl.cols()));
    }

    ML_CVLOG(5, this->get_logger_id(), "QI:");
    for (int row = 0; row < this->_q_delta_impl.rows(); ++row) {
      ML_CVLOG(5, this->get_logger_id(),
               "  " << this->_q_delta_impl.block(row, 0, 1, this->_q_delta_impl.cols()));
    }
  }
}  // ::pfasst
