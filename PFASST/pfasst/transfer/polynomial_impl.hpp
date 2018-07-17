#include "pfasst/transfer/polynomial.hpp"

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <vector>
using std::shared_ptr;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"
#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate_initial(const shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                                                                   shared_ptr<typename TransferTraits::fine_sweeper_t> fine)
  {
    ML_CVLOG(1, "TRANS", "interpolate initial value only");

    auto& coarse_factory = coarse->get_encap_factory();
    auto& fine_factory = fine->get_encap_factory();

    // c_delta = restrict(u_0^F) - u_0^C
    auto coarse_delta = coarse_factory.create();
    // c_delta = restrict(u_0^F)
    this->restrict_data(fine->get_initial_state(), coarse_delta);
    // c_delta -= c_0
    coarse_delta->scaled_add(-1.0, coarse->get_initial_state());

    // f_delta = interpolate(c_delta)
    auto fine_delta = fine_factory.create();
    this->interpolate_data(coarse_delta, fine_delta);

    // u_0^F = u_0^F - f_delta
    fine->initial_state()->scaled_add(-1.0, fine_delta);

    // update function evaluations on fine level's initial state
    fine->reevaluate(true);
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate(const shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                                                           shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                                           const bool initial)
  {
    ML_CVLOG(1, "TRANS", "interpolate");

    if (coarse->get_quadrature()->left_is_node() || fine->get_quadrature()->left_is_node()) {
      ML_CLOG(ERROR, "TRANS", "Time interpolation with left time point as a node is still not supported.");
      throw std::runtime_error("time interpolation with left time point as node");
    }

    ML_CLOG_IF(fine->get_quadrature()->get_num_nodes() != coarse->get_quadrature()->get_num_nodes(),
            WARNING, "TRANS", "interpolation between different number of nodes not tested!");

    if (initial) {
      this->interpolate_initial(coarse, fine);
    }

    this->setup_tmat(fine->get_quadrature(), coarse->get_quadrature());
//     ML_CVLOG(1, "TRANS", "tmat: " << this->tmat);

    // +1 here for additional value in states
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes() + 1;

    auto& coarse_factory = coarse->get_encap_factory();
    auto& fine_factory = fine->get_encap_factory();

    // step 1: compute coarse level correction
    std::vector<shared_ptr<typename traits::fine_encap_t>> fine_deltas(num_coarse_nodes);
    std::generate(fine_deltas.begin(), fine_deltas.end(),
             [&fine_factory]() { return fine_factory.create(); });

    //  c_delta = restrict(u_m^F) - u_m^C
    //  f_delta = interpolate(c_delta)
    auto coarse_delta = coarse_factory.create();
    for (size_t m = 1; m < num_coarse_nodes; ++m) {
      this->restrict_data(fine->get_states()[m], coarse_delta);
      coarse_delta->scaled_add(-1.0, coarse->get_states()[m]);
//       ML_CVLOG(1, "TRANS", "  cd["<<m<<"]: " << to_string(coarse_delta));
      this->interpolate_data(coarse_delta, fine_deltas[m]);
//       ML_CVLOG(1, "TRANS", "  fd["<<m<<"]: " << to_string(fine_deltas[m]));
    }

    // step 2: add coarse level correction onto fine level's states
//     ML_CVLOG(1, "TRANS", "fine states and deltas before interpolation:");
//     for (size_t m = 0; m < num_coarse_nodes; ++m) {
//       ML_CVLOG(1, "TRANS", "  f["<<m<<"]:  " << to_string(fine->get_states()[m]));
//       ML_CVLOG(1, "TRANS", "  fd["<<m<<"]: " << to_string(fine_deltas[m]));
//     }

    // u^F = u^F - T f_delta
    encap::mat_apply(fine->states(), -1.0, this->tmat, fine_deltas, false);

//     ML_CVLOG(1, "TRANS", "fine states after interpolation:");
//     for (auto& n : fine->states()) {
//       ML_CVLOG(1, "TRANS", "  " << to_string(n));
//     }

    // update function evaluations on fine level
    fine->reevaluate();
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                                                shared_ptr<typename TransferTraits::fine_encap_t> fine)
  {
    UNUSED(coarse); UNUSED(fine);
    throw std::runtime_error("interpolation for generic Encapsulations");
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict_initial(const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                                                shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse)
  {
    ML_CVLOG(1, "TRANS", "restrict initial value only");
    // M * fine->get_initial_state()
    shared_ptr<typename TransferTraits::fine_encap_t> M_initial_state= fine->get_encap_factory().create();
    fine->get_initial_state()->apply_Mass(*(fine->get_M_dune()) , M_initial_state);
    this->restrict_u(M_initial_state , coarse->_M_initial);    
    
    this->restrict_data(fine->get_initial_state(), coarse->initial_state());


  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict(const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                                        shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                                                        const bool initial)
  {
    ML_CVLOG(1, "TRANS", "restrict");

    if (coarse->get_quadrature()->left_is_node() || fine->get_quadrature()->left_is_node()) {
      ML_CLOG(ERROR, "TRANS", "Time restriction with left time point as a node is still not supported.");
      throw std::runtime_error("time restriction with left time point as node");
    }


    if (initial) {
      this->restrict_initial(fine, coarse);
    }

    // +1 here for additional value in states
    const auto coarse_nodes = coarse->get_quadrature()->get_nodes();
    const auto fine_nodes = fine->get_quadrature()->get_nodes();
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes() + 1;

    // this commented out stuff is probably required for non-equal sets of time nodes
//     const int factor = ((int)num_fine_nodes - 1) / ((int)num_coarse_nodes - 1);

    for (size_t m = 1; m < num_coarse_nodes; ++m) {
//       if (coarse_nodes[m] != fine_nodes[m * factor]) {
//         CLOG(ERROR, "TRANS") << "coarse nodes are not nested within fine ones."
//                              << "coarse: " << coarse_nodes << " fine: " << fine_nodes;
//         throw NotImplementedYet("non-nested nodes");
//       }
      this->restrict_data(fine->get_states()[m], coarse->states()[m]);
    }

    coarse->reevaluate();
  }

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                                             shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
  {
    UNUSED(coarse); UNUSED(fine);
    throw std::runtime_error("restriction for generic Encapsulations");
  }
  
  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::restrict_u(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                                             shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
  {
    UNUSED(coarse); UNUSED(fine);
    throw std::runtime_error("restriction for generic Encapsulations");
  }
  

  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::fas(const typename TransferTraits::fine_time_t& dt,
                                                   const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                                   shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse)
  {
    ML_CVLOG(1, "TRANS", "compute FAS correction");

    const auto coarse_nodes = coarse->get_quadrature()->get_nodes();
    const auto fine_nodes = fine->get_quadrature()->get_nodes();
    const size_t num_coarse_nodes = coarse->get_quadrature()->get_num_nodes();
    const size_t num_fine_nodes = fine->get_quadrature()->get_num_nodes();

    if (num_fine_nodes != num_coarse_nodes) {
      ML_CLOG(ERROR, "TRANS", "FAS Correction for different time scales not supported yet.");
      throw std::runtime_error("fas correction for different time scales");
    } else if (!equal(coarse_nodes.cbegin(), coarse_nodes.cend(), fine_nodes.cbegin())) {
      ML_CLOG(ERROR, "TRANS", "FAS Correction for different time nodes not supported yet.");
      throw std::runtime_error("fas correction for different sets of nodes");
    }
    assert(coarse_nodes.size() == fine_nodes.size());

    auto& coarse_factory = coarse->get_encap_factory();

    std::vector<shared_ptr<typename traits::coarse_encap_t>> fas(num_coarse_nodes + 1);
    std::generate(fas.begin(), fas.end(), [&coarse_factory]() { return coarse_factory.create(); });


    //std::cout << "coasre integrate " << num_coarse_nodes<< std::endl;
    
    const auto coarse_integral = coarse->integrate(dt);  

    /*for (size_t m = 1; m < num_coarse_nodes; ++m) {
      shared_ptr<typename TransferTraits::coarse_encap_t> coarse_u= coarse->get_encap_factory().create();	
      this->restrict_data(fine->get_states()[m], coarse_u);  
      //coarse->_impl_rhs_restrict[m] =  coarse->evaluate_rhs_impl(0, coarse_u);
      
      //for(int i=0; i<coarse->_impl_rhs_restrict[m]->data().size(); i++){
        //std::cout << " r rt " << coarse->_impl_rhs_restrict[m]->data()[i] << " " << coarse->_impl_rhs[m]->data()[i] << std::endl;
      //}  
    }*/
    
    for (size_t m = 0; m < num_coarse_nodes + 1; ++m) {
      
      shared_ptr<typename TransferTraits::coarse_encap_t> coarse_u= coarse->get_encap_factory().create();	
      this->restrict_data(fine->get_states()[m], coarse_u);  

      //coarse->get_M_dune()->mv(coarse_u->get_data(),  coarse_integral[m]->data());
      coarse_u->apply_Mass(*(coarse->get_M_dune()), coarse_integral[m]);

      coarse_integral[m]->data() *= -1;
      
      coarse_integral[m]->scaled_add(1,  coarse->integrate(dt)[m]);
    }
 
    //std::cout << "fine integrate " << num_fine_nodes << std::endl;
 
    const auto fine_integral = fine->integrate(dt);
    
    
    
    //this->restrict_u(fine->get_M_dune(), ); //raus
    
    for (size_t m = 0; m < num_fine_nodes + 1; ++m) {
      //fine->get_M_dune()->mv(fine->get_states()[m]->get_data(),  fine_integral[m]->data());
      fine->get_states()[m]->apply_Mass(*(fine->get_M_dune()),  fine_integral[m]);
      fine_integral[m]->data() *= -1;
      
      
      fine_integral[m]->scaled_add(1,  fine->integrate(dt)[m]);
    } 
    


    for (size_t m = 0; m < num_coarse_nodes + 1; ++m) {
      this->restrict_u(fine_integral[m], fas[m]);
      fas[m]->scaled_add(-1.0, coarse_integral[m]);
      coarse->tau()[m]->data() = fas[m]->get_data();
      //coarse->tau()[m]->data() *=0;
      /*for(int i=0; i<coarse->tau()[m]->data().size(); i++){
      std::cout << coarse->tau()[m]->data() << std::endl;
      }*/
      //std::cout << "norm tau" << coarse->tau()[m]->norm0() << std::endl;
      //coarse->tau()[m]->data() *=-1;
    }
  }


  template<class TransferTraits, typename Enabled>
  void
  PolynomialTransfer<TransferTraits, Enabled>::setup_tmat(const shared_ptr<quadrature::IQuadrature<typename TransferTraits::fine_time_t>> fine_quad,
                                                          const shared_ptr<quadrature::IQuadrature<typename TransferTraits::coarse_time_t>> coarse_quad)
  {
    if (this->tmat.rows() == 0) {
      auto coarse_nodes = coarse_quad->get_nodes();
      auto fine_nodes = fine_quad->get_nodes();

      coarse_nodes.insert(coarse_nodes.begin(), 0.0);
      fine_nodes.insert(fine_nodes.begin(), 0.0);

      ML_CLOG_IF(!equal(coarse_nodes.cbegin(), coarse_nodes.cend(), fine_nodes.cbegin()),
              WARNING, "TRANS", "interpolation for different sets of nodes not tested.");

      this->tmat = quadrature::compute_interp<typename TransferTraits::fine_time_t>(coarse_nodes,
                                                                                       fine_nodes);
    }
  }
}  // ::pfasst
