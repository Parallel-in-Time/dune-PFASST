#include "pfasst/controller/two_level_mlsdc_n.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <type_traits>
using std::shared_ptr;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT, class CommT>
  TwoLevelMLSDC<TransferT, CommT>::TwoLevelMLSDC()
    : Controller<TransferT, CommT>()
  {
    TwoLevelMLSDC<TransferT, CommT>::init_loggers();
    this->set_logger_id("MLSDC");
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::init_loggers()
  {
    log::add_custom_logger("MLSDC");
    log::add_custom_logger("LVL_COARSE");
    log::add_custom_logger("LVL_FINE");
  }

  template<class TransferT, class CommT>
  size_t
  TwoLevelMLSDC<TransferT, CommT>::get_num_levels() const
  {
    size_t num = 0;
    if (this->_coarse_level != nullptr) {
      num++;
    }
    if (this->_fine_level != nullptr) {
      num++;
    }
    return num;
  }

  template<class TransferT, class CommT>
  template<class SweeperT>
  void
  TwoLevelMLSDC<TransferT, CommT>::add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse)
  {

      if (std::is_same<SweeperT, typename transfer_t::traits::coarse_sweeper_t>::value) {
        this->_coarse_level = sweeper;
        this->get_coarse()->set_logger_id("LVL_COARSE");
      } else {
        ML_CLOG(ERROR, this->get_logger_id(), "Type of given Sweeper ("
          << typeid(SweeperT).name() << ") is not applicable as Coarse Sweeper ("
          << typeid(typename transfer_t::traits::coarse_sweeper_t).name() << ").");
        throw std::logic_error("given sweeper can not be used as coarse sweeper");
      }
    
    
  }
  
  
  template<class TransferT, class CommT>
  template<class SweeperT>
  void
  TwoLevelMLSDC<TransferT, CommT>::add_sweeper(shared_ptr<SweeperT> sweeper)
  {
   

      if (std::is_same<SweeperT, typename transfer_t::traits::fine_sweeper_t>::value) {
        
        
        this->_fine_level = sweeper;
        this->get_fine()->set_logger_id("LVL_FINE");
        
        
        
      } else {
        ML_CLOG(ERROR, this->get_logger_id(), "Type of given Sweeper ("
          << typeid(SweeperT).name() << ") is not applicable as Fine Sweeper ("
          << typeid(typename transfer_t::traits::fine_sweeper_t).name() << ").");
        throw std::logic_error("given sweeper can not be used as fine sweeper");
      }
    
  }
  
  

  template<class TransferT, class CommT>
  const shared_ptr<typename TransferT::traits::coarse_sweeper_t>
  TwoLevelMLSDC<TransferT, CommT>::get_coarse() const
  {
    return this->_coarse_level;
  }

  template<class TransferT, class CommT>
  shared_ptr<typename TransferT::traits::coarse_sweeper_t>
  TwoLevelMLSDC<TransferT, CommT>::get_coarse()
  {
    return this->_coarse_level;
  }

  template<class TransferT, class CommT>
  const shared_ptr<typename TransferT::traits::fine_sweeper_t>
  TwoLevelMLSDC<TransferT, CommT>::get_fine() const
  {
    return this->_fine_level;
  }

  template<class TransferT, class CommT>
  shared_ptr<typename TransferT::traits::fine_sweeper_t>
  TwoLevelMLSDC<TransferT, CommT>::get_fine()
  {
    return this->_fine_level;
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::set_options()
  {
    Controller<TransferT, CommT>::set_options();

    this->get_fine()->set_options();
    this->get_coarse()->set_options();
  }


  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::setup()
  {
    assert(this->get_transfer() != nullptr);

    Controller<TransferT, CommT>::setup();

    if (this->get_num_levels() != 2) {
      ML_CLOG(ERROR, this->get_logger_id(), "Two levels (Sweeper) must have been added for Two-Level-MLSDC.");
      throw std::logic_error("Two-Level-MLSDC requires two levels");
    }

    ML_CVLOG(1, this->get_logger_id(), "setting up coarse level");
    this->get_coarse()->status() = this->get_status();
    this->get_coarse()->setup();

    ML_CVLOG(1, this->get_logger_id(), "setting up fine level");
    this->get_fine()->status() = this->get_status();
    this->get_fine()->setup();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::run()
  {
    Controller<TransferT, CommT>::run();

    do {
      ML_CLOG(INFO, this->get_logger_id(), "");
      ML_CLOG(INFO, this->get_logger_id(), "Time Step " << (this->get_status()->get_step() + 1)
                                        << " of " << this->get_status()->get_num_steps());

      this->status()->set_primary_state(PrimaryState::PREDICTING);

      bool last_predict=false;
      // iterate on each time step
      do {
        if (this->get_status()->get_primary_state() == (+PrimaryState::PREDICTING)) {
          ML_CLOG(INFO, this->get_logger_id(), "");
          ML_CLOG(INFO, this->get_logger_id(), "Iteration 0 (MLSDC Prediction)");

          assert(this->get_status()->get_iteration() == 0);


          //this->get_transfer()->restrict_initial(this->get_fine(), this->get_coarse());



          //this->get_coarse()->spread_Newton(); 
	  //this->get_fine()->spread_Newton();
          //this->get_coarse()->save(); this->get_fine()->save();

          this->get_fine()->predict();

          this->cycle_down();
          this->sweep_coarse();
          this->cycle_up();
          this->sweep_fine();


        } else {
          ML_CLOG(INFO, this->get_logger_id(), "");
          ML_CLOG(INFO, this->get_logger_id(), "Iteration " << this->get_status()->get_iteration());


          this->cycle_down();
          this->sweep_coarse();
          this->cycle_up();
          this->sweep_fine();

        }

        this->status()->set_primary_state(PrimaryState::INTER_ITER);
      } while(this->advance_iteration());

      /*auto grid = (*(this->get_fine())).get_grid();
      typedef Dune::BlockVector<Dune::FieldVector<double, 2> > VectorType;
      typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ColumnType;
      typedef Dune::YaspGrid<2> GridType; //ruth_dim
      typedef GridType::LevelGridView GridView;
      GridType::LevelGridView gridView = grid->levelGridView(1);
      Dune::VTKWriter<GridView> vtkWriter(gridView);


      VectorType x = this->get_fine()->get_end_state()->data();


      stringstream toss;
      toss << this->get_status()->get_step();
      string name = toss.str();
      //string name = boost::lexical_cast<string>(this->get_status()->get_step());

      ColumnType sol_u, sol_v;
      sol_u.resize(x.size());
      sol_v.resize(x.size());
      for (int i =0; i< x.size(); ++i)
      {
        sol_u[i] = x[i][0];
        sol_v[i] = x[i][1];
      }


      vtkWriter.addVertexData(sol_u, "fe_solution_u");
      vtkWriter.addVertexData(sol_v, "fe_solution_v");


      vtkWriter.write("gray_scott" + name);*/



    } while(this->advance_time());
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelMLSDC<TransferT, CommT>::advance_time(const size_t& num_steps)
  {
    if (Controller<TransferT, CommT>::advance_time(num_steps)) {
      this->get_fine()->advance(num_steps);
      this->get_coarse()->advance(num_steps);
      return true;
    } else {
      return false;
    }
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelMLSDC<TransferT, CommT>::advance_time()
  {
    return this->advance_time(1);
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelMLSDC<TransferT, CommT>::advance_iteration() //hier jetzt mit alternativem residuumscheck...
  {
    this->status()->set_secondary_state(SecondaryState::CONV_CHECK);
    std::cout << "check if converged" << std::endl;
    this->get_coarse()->converged(false);
    std::cout << "check if converged1" << std::endl;
    if (this->get_fine()->converged(false)) {
      ML_CLOG(INFO, this->get_logger_id(), "FINE sweeper has converged.");
      this->status()->set_primary_state(PrimaryState::CONVERGED);
      return false;

    } else if (Controller<TransferT, CommT>::advance_iteration()) {
      ML_CLOG(INFO, this->get_logger_id(), "FINE sweeper has not yet converged and additional iterations to do.");
      this->get_fine()->save();
      this->get_coarse()->save();
      this->status()->set_primary_state(PrimaryState::ITERATING);
      return true;

    } else {
      ML_CLOG(INFO, this->get_logger_id(), "FINE sweeper has not yet converged and no more iterations to do.");
      this->status()->set_primary_state(PrimaryState::FAILED);
      return false;
    }
        std::cout << "check if converged done" << std::endl;
  }


  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::predict_coarse()
  {
    ML_CLOG(INFO, this->get_logger_id(), "Predicting on COARSE level");

    this->status()->set_secondary_state(SecondaryState::PRE_ITER_COARSE);
    this->get_coarse()->pre_predict();

    this->status()->set_secondary_state(SecondaryState::ITER_COARSE);
    this->get_coarse()->predict();

    this->status()->set_secondary_state(SecondaryState::POST_ITER_COARSE);
    this->get_coarse()->post_predict();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::predict_fine()
  {
    ML_CLOG(INFO, this->get_logger_id(), "Predicting on FINE level");

    this->status()->set_secondary_state(SecondaryState::PRE_ITER_FINE);
    this->get_fine()->pre_predict();

    this->status()->set_secondary_state(SecondaryState::ITER_FINE);
    this->get_fine()->predict();

    this->status()->set_secondary_state(SecondaryState::POST_ITER_FINE);
    this->get_fine()->post_predict();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::sweep_coarse()
  {
    ML_CLOG(INFO, this->get_logger_id(), "Sweeping on COARSE level");

    this->status()->set_secondary_state(SecondaryState::PRE_ITER_COARSE);
    this->get_coarse()->pre_sweep();

    this->status()->set_secondary_state(SecondaryState::ITER_COARSE);
    this->get_coarse()->sweep();

    this->status()->set_secondary_state(SecondaryState::POST_ITER_COARSE);
    this->get_coarse()->post_sweep();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::sweep_fine()
  {
    ML_CLOG(INFO, this->get_logger_id(), "Sweeping on FINE level");

    this->status()->set_secondary_state(SecondaryState::PRE_ITER_FINE);
    this->get_fine()->pre_sweep();

    this->status()->set_secondary_state(SecondaryState::ITER_FINE);
    this->get_fine()->sweep();

    this->status()->set_secondary_state(SecondaryState::POST_ITER_FINE);
    this->get_fine()->post_sweep();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::cycle_down()
  {
    ML_CLOG(INFO, this->get_logger_id(), "Restrict onto coarse level");
    this->status()->set_secondary_state(SecondaryState::CYCLE_DOWN);

    this->get_transfer()->restrict(this->get_fine(), this->get_coarse(), true);
    this->get_transfer()->fas(this->get_status()->get_dt(), this->get_fine(), this->get_coarse());
    this->get_coarse()->save();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelMLSDC<TransferT, CommT>::cycle_up()
  {
    ML_CLOG(INFO, this->get_logger_id(), "Interpolate onto fine level");
    this->status()->set_secondary_state(SecondaryState::CYCLE_UP);
    std::cout << "666666666666666666666666666666666666666666666666666666666666666666666666666 vor interpolate" << std::endl;	
    this->get_transfer()->interpolate(this->get_coarse(), this->get_fine(), true);
    std::cout << "666666666666666666666666666666666666666666666666666666666666666666666666666 nach interpolate" << std::endl;	
  }
}  // ::pfasst
