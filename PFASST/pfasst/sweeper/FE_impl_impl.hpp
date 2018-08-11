
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

    const auto num_nodes = this->get_quadrature()->get_num_nodes();
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

    for (size_t m = 0; m < num_nodes; ++m) {
      this->states()[m + 1]->data() = this->states()[m]->data();	
      this->_impl_rhs[m + 1] = this->evaluate_rhs_impl(tm, this->get_states()[m + 1]);
      tm += dt *  (nodes[m+1] - nodes[m]);
      ML_CVLOG(1, this->get_logger_id(), "");

    }       	

  }

  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::post_predict()
  {
    Sweeper<SweeperTrait, BaseFunction, Enabled>::post_predict();
  }

  
  
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::pre_sweep()
  {
 
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
    ML_CVLOG(6, this->get_logger_id(), "           += dt * Q * f_im");

    encap::mat_apply(this->_q_integrals, dt, q_mat, this->_impl_rhs, true); //false

    ML_CVLOG(4, this->get_logger_id(), "  subtracting function evaluations of previous iteration and adding FAS correction");


    for (size_t m = 0; m < num_nodes; ++m) {
      for (size_t n = 0; n < m + 1; ++n) {

        this->_q_integrals[m + 1]->scaled_add(-dt * this->_q_delta_impl(m + 1, n + 1), this->_impl_rhs[n + 1]);
      }

      this->_q_integrals[m + 1]->scaled_add(1.0, this->get_tau()[m + 1]);
    }

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

    typename traits::time_t tm = t;
    for (size_t m = 0; m < num_nodes; ++m) {
      ML_CVLOG(4, this->get_logger_id(), "propagating from t["<<m<<"]=" << (t + (dt * nodes[m]))
                                                   << " to t["<<(m+1)<<"]=" << (t + (dt * nodes[m+1])));
	
    shared_ptr<typename traits::encap_t> rhs = this->get_encap_factory().create();

    if (is_coarse){
        rhs->data() =  this->_M_initial->get_data(); 
    }else{
	this->get_states().front()->apply_Mass(M_dune, rhs);  
    }
    // rhs += dt * \sum_{i=0}^m (QI_{m+1,i} fI(u_i^{k+1}) + QE_{m+1,i-1} fE(u_{i-1}^{k+1}) ) + QE_{m+1,m} fE(u_{m}^{k+1})
    for (size_t n = 0; n <= m; ++n) {
        rhs->scaled_add(dt * this->_q_delta_impl(m + 1, n), this->_impl_rhs[n]);
    } 

    rhs->scaled_add(1.0, this->_q_integrals[m + 1]);

    // solve the implicit part
    ML_CVLOG(4, this->get_logger_id(), "  solve(u["<<(m+1)<<"] - dt * QI_{"<<(m+1)<<","<<(m+1)<<"} * f_im["<<(m+1)<<"] = rhs)");
    this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, dt * this->_q_delta_impl(m+1, m+1), rhs);

    tm += dt * this->_q_delta_impl(m+1, m+1);
      
    ML_CVLOG(4, this->get_logger_id(), "");


    }

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
      assert(this->_impl_rhs.front() != nullptr && this->_impl_rhs.back() != nullptr);
      this->_impl_rhs.front()->data() = this->_impl_rhs.back()->get_data();
    } else {
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
      this->_impl_rhs.front() = this->evaluate_rhs_impl(t0, this->get_initial_state());

    } else {
      const typename traits::time_t dt = this->get_status()->get_dt();
      auto nodes = this->get_quadrature()->get_nodes();
      nodes.insert(nodes.begin(), 0.0);
      for (size_t m = 0; m < this->get_quadrature()->get_num_nodes() + 1; ++m) {
        const typename traits::time_t t = t0 + dt * nodes[m];
        assert( this->_impl_rhs[m] != nullptr);
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
    auto result = encap::mat_mul_vec(dt, q_mat, this->_impl_rhs);
    return result;
  }
  
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_t>>
  IMEX<SweeperTrait, BaseFunction, Enabled>::integrate_new(const typename SweeperTrait::time_t& dt)
  {
    auto const q_mat = this->get_quadrature()->get_q_mat();
    auto result = encap::mat_mul_vec(dt, q_mat, this->_impl_rhs_restrict);
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

      this->end_state()->data() = this->get_initial_state()->get_data();
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_impl_rhs)[0]);
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
    const typename traits::time_t dt = this->get_status()->get_dt();
    const size_t num_nodes = this->get_quadrature()->get_num_nodes() + 1;
    if (only_last) {
      const size_t cols = this->get_quadrature()->get_q_mat().cols();
      const size_t rows = this->get_quadrature()->get_q_mat().rows();
          
      if (is_coarse){

        this->residuals().back()->data() =  this->_M_initial->get_data(); 
      }else{
	this->get_initial_state()->apply_Mass(M_dune, this->residuals().back());
      }
      shared_ptr<typename traits::encap_t> uM = this->get_encap_factory().create();	
      this->get_states().back()->apply_Mass(M_dune,uM);	
      this->get_states().back()->apply_Mass(M_dune, uM);	
      this->residuals().back()->scaled_add(-1.0, uM);	

      this->residuals().back()->scaled_add(1.0, this->get_tau().back());
      for (size_t n = 0; n < cols; ++n) {
        this->residuals().back()->scaled_add(dt * this->get_quadrature()->get_q_mat()(rows - 1, n), this->_impl_rhs[n]);
      }
    } else {

      for (size_t m = 0; m < num_nodes; ++m) {
        assert(this->get_states()[m] != nullptr);
        assert(this->residuals()[m] != nullptr);
    
    	if (is_coarse){
        	this->residuals()[m]->data() =  this->_M_initial->get_data(); 
    	}else{
		this->get_initial_state()->apply_Mass(M_dune,  this->residuals()[m]);
    	}   
    
	shared_ptr<typename traits::encap_t> uM = this->get_encap_factory().create();
	this->get_states()[m]->apply_Mass(M_dune, uM);
	this->residuals()[m]->scaled_add(-1.0,uM);
        assert(this->get_tau()[m] != nullptr);
        
        this->residuals()[m]->scaled_add(1.0, this->get_tau()[m]);
        
      }

      ML_CVLOG(5, this->get_logger_id(), "vor  res += dt * Q * F_im");

      encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_impl_rhs, false); //DAS IST FALSCH RUTH!!! 
      	    
    }
  }


  
  
  template<class SweeperTrait, class BaseFunction, typename Enabled>
  void
  IMEX<SweeperTrait, BaseFunction, Enabled>::compute_residuals()
  {
    this->compute_residuals(false);
  }



  template<class SweeperTrait, class BaseFunction, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_t>
  IMEX<SweeperTrait, BaseFunction, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u)
  {
    UNUSED(t); UNUSED(u);
    throw std::runtime_error("evaluation of implicit part of right-hand-side");
  }

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


    }else{



        for (size_t m = 1; m < num_nodes + 1; ++m) {
            for (size_t n = m; n < num_nodes + 1; ++n) {
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
