#include "dune_vec.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/logging.hpp"
#include "pfasst/util.hpp"


namespace pfasst
{
  namespace encap
  {
    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(const size_t size)
      : _data(size) , size(size)
    {
      this->zero();
    }

    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(const typename EncapsulationTrait::data_t& data)
      : Encapsulation<EncapsulationTrait>(data.size())
    {
      this->data() = data;
    }

    template<class EncapsulationTrait>
    Encapsulation<EncapsulationTrait>&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::operator=(const typename EncapsulationTrait::data_t& data)
    {
      this->data() = data;
      return *this;
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::data_t&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::data()
    {
      return this->_data;
    }

    template<class EncapsulationTrait>
    const typename EncapsulationTrait::data_t&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_data() const
    {
      return this->_data;
    }

    template<class EncapsulationTrait>
    size_t
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_total_num_dofs() const
    {
      return this->get_data().size();
    }

    /*template<class EncapsulationTrait>
    std::array<int, EncapsulationTrait::DIM>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_dimwise_num_dofs() const
    {
      std::array<int, EncapsulationTrait::DIM> dimwise_ndofs;
      switch (EncapsulationTrait::DIM) {
        case 1:
          dimwise_ndofs.fill((int)this->get_total_num_dofs()); //setzt 1,2 v 3 dimensionen und quadratischen raum vorraus
          break;
        case 2:
          dimwise_ndofs.fill((int)sqrt(this->get_total_num_dofs()));
          break;
        case 3:
          dimwise_ndofs.fill((int)cbrt(this->get_total_num_dofs()));
          break;
        default:
          ML_CLOG(FATAL, "ENCAP", "unsupported spatial dimension: " << EncapsulationTrait::DIM);
          throw std::runtime_error("unsupported spatial dimension");
      }

      return dimwise_ndofs;
    }*/

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::zero()
    {
      std::fill(this->data().begin(), this->data().end(), typename EncapsulationTrait::spatial_t(0.0));
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::scaled_add(const typename EncapsulationTrait::time_t& a,
                                   const shared_ptr<Encapsulation<EncapsulationTrait>> y)
    {

      assert(this->get_data().size() == y->data().size());

      std::transform(this->get_data().begin(), this->get_data().end(), y->data().begin(),
                this->data().begin(),
                [a](const typename EncapsulationTrait::spatial_t& xi,
                    const typename EncapsulationTrait::spatial_t& yi) { return xi + a * yi; });

    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::spatial_t
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::norm0() const
    {
    
      /*int my_rank, num_pro;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
      MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
      MPI_Comm comm_x, comm_t; 
      int myid, xcolor, tcolor;
      int space_num=2;
      xcolor = (my_rank / space_num);
      tcolor = my_rank % space_num;
      std::cout << "------------------------------------------------------------------------------ vor split im norm() " << my_rank << std::endl;
      //MPI_Comm_split( MPI_COMM_WORLD, xcolor, my_rank, &comm_x );
      MPI_Comm_split( MPI_COMM_WORLD, tcolor, my_rank, &comm_t );
      std::cout << "------------------------------------------------------------------------------nach split im norm() " << my_rank << std::endl;
      //std::cout << " in der normberechnung " << std::endl;*/
      double max = std::abs(*(std::max_element(this->get_data().begin(), this->get_data().end(),
                               [](const typename EncapsulationTrait::spatial_t& a,
                                  const typename EncapsulationTrait::spatial_t& b)
                                 { return std::abs(a) < std::abs(b); })));
      double global_max=max;
      //std::cout << "------------------------------------------------------------------------------ neue norm" << max << std::endl;
      //std::exit(0);
      //MPI_Allreduce(&max,&global_max,1,MPI_DOUBLE,MPI_MAX,comm_t);
      //std::cout << "------------------------------------------------------------------------------ neue norm all" << max << std::endl;
      return global_max;
      /*using Norm =  EnergyNorm<MatrixType,VectorType>;
      auto parallel_energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(this->A_dune, restrictToMaster, gridView.grid().comm());
      return parallel_energyNorm(this->get_data());*/
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::spatial_t
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::norm0(bool ignore, MPI_Comm comm) const
    {
    
      //return this->norm0();
      std::cout << "in norm0 mpi";
      /*int my_rank, num_pro;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
      MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
      MPI_Comm comm_x, comm_t; 
      int myid, xcolor, tcolor;
      int space_num=2;
      xcolor = (my_rank / space_num);
      tcolor = my_rank % space_num;*/
      int rank, num_pro;
      MPI_Comm_rank(comm, &rank );
      MPI_Comm_size(comm, &num_pro );
      //std::cout << "------------------------------------------------------------------------------ vor split im norm() " << rank << std::endl;
      //MPI_Comm_split( MPI_COMM_WORLD, xcolor, my_rank, &comm_x );
      //MPI_Comm_split( MPI_COMM_WORLD, tcolor, my_rank, &comm_t );
      //std::cout << "------------------------------------------------------------------------------nach split im norm() " << my_rank << std::endl;
      
      
      double max;
      auto begin=this->get_data().begin();begin++; begin++;begin++;
      auto end=this->get_data().end();end--; end--;end--;
      if(rank==0) {begin--;begin--;begin--;}
      if(rank==num_pro-1){end++;end++;end++;} 
      max = std::abs(*(std::max_element(begin, end,
                               [](const typename EncapsulationTrait::spatial_t& a,
                                  const typename EncapsulationTrait::spatial_t& b)
                                 { return std::abs(a) < std::abs(b); })));
      double global_max=max;
      //std::cout << "------------------------------------------------------------------------------ neue norm ausgerechnet " << max << " " << rank <<std::endl;
      MPI_Allreduce(&max,&global_max,1,MPI_DOUBLE,MPI_MAX,comm);
      return global_max;
      
      /*using Norm =  EnergyNorm<MatrixType,VectorType>;
      auto parallel_energyNorm = Dune::ParMG::parallelEnergyNorm<VectorType>(this->A_dune, restrictToMaster, gridView.grid().comm());
      return parallel_energyNorm(this->get_data());*/
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::apply_Mass(typename traits::mass_t mass, shared_ptr<Encapsulation<EncapsulationTrait>> sol)//, EncapsulationTrait &sol)
    {
    	//typename EncapsulationTrait::data_t copy= this->get_data();
	mass.mv(this->data(), sol->data());//sol->data());  
    } 

    template<class EncapsulationTrait>
    template<class CommT>
    bool
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::probe(shared_ptr<CommT> comm, const int src_rank, const int tag)
    {
      return comm->probe(src_rank, tag);
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::send(shared_ptr<CommT> comm, const int dest_rank,
                              const int tag, const bool blocking)
    {
      ML_CVLOG(2, "ENCAP", "sending data: " << this->get_data());
      if (blocking) {
        //comm->send(this->get_data().data(), this->get_data().size(), dest_rank, tag);
        comm->send(&(this->data()[0][0]), this->size, dest_rank, tag);
      } else {
        //comm->isend(this->get_data().data(), this->get_data().size(), dest_rank, tag);
        comm->isend(&(this->data()[0][0]), this->size, dest_rank, tag);
      }
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::recv(shared_ptr<CommT> comm, const int src_rank,
                              const int tag, const bool blocking)
    {
      if (blocking) {
        //comm->recv(this->data().data(), this->get_data().size(), src_rank, tag);
        comm->recv(&(this->data()[0][0]), this->size, src_rank, tag);
      } else {
        //comm->irecv(this->data().data(),    this->get_data().size(), src_rank, tag);
        comm->irecv(&(this->data()[0][0]), this->size, src_rank, tag);
      }
      ML_CVLOG(2, "ENCAP", "received data: " << this->get_data());
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::bcast(shared_ptr<CommT> comm, const int root_rank)
    {
      //comm->bcast(this->data().data(), this->get_data().size(), root_rank);
      comm->bcast(&(this->data()[0][0]), this->size, root_rank);
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::log(el::base::type::ostream_t& os) const
    {
      os << "not implemented";// "FieldVector" << pfasst::join(this->get_data(), ", ");
    }


    template<class EncapsulationTrait>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::EncapsulationFactory(const size_t size)
      : _size(size)
    {}

    /*template<class EncapsulationTrait>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::EncapsulationFactory(typename EncapsulationTrait::gfs_t gfs)
      : _gfs(gfs)
    {}*/

    template<class EncapsulationTrait>
    shared_ptr<Encapsulation<EncapsulationTrait>>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::create() const
    {
      return std::make_shared<Encapsulation<EncapsulationTrait>>(this->size());
    }

    /*template<class EncapsulationTrait>
    void
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::set_gfs(typename EncapsulationTrait::gfs_t& gfs)
    {
      this->_gfs = std::make_shared<typename EncapsulationTrait::gfs_t>(gfs);
    }*/

    template<class EncapsulationTrait>
    void
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::set_size(const size_t& size)
    {
      this->_size = size;
    }
    
            
    /*template<class EncapsulationTrait>
    void
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::set_FE_manager(std::shared_ptr<fe_manager> FinEl, int nlevel, MatrixType A_dune)
    {
    	int rank, num_pro;
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank );
    	MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        auto basis = FinEl->get_basis(nlevel);
        auto grid = FinEl->get_grid();
        auto bs = basis->size();
        
  	using MGSetup = Dune::ParMG::ParallelMultiGridSetup< BasisFunction, MatrixType, VectorType >;
  	MGSetup mgSetup{*grid,grid->maxLevel() - (nlevel)};
  	auto gridView = mgSetup.bases_.back().gridView();
 	using MG = Dune::ParMG::Multigrid<VectorType>;
  	MG mg;
  
    	using namespace Dune::ParMG;
    	auto& levelOp = mgSetup.levelOps_;
    	auto df_pointer = std::make_shared<MatrixType>(A_dune);	
    	mgSetup.matrix(df_pointer);
    	auto fineIgnore = std::make_shared< Dune::BitSetVector<1> >(bs);
    	for (std::size_t i = 0; i < bs; ++i){
      		(*fineIgnore)[i] = false;
      		if(i==0 &&rank==0) (*fineIgnore)[i] = true;
      		if(rank==num_pro - 1 && i== bs-1) (*fineIgnore)[i] = true;
      	}
      	    		//std::cout << "nach ignore gesetzt " << std::endl;

    	mgSetup.ignore(fineIgnore);
    	mgSetup.setupLevelOps();
    	double dampening =1.0;
    	mgSetup.setupSmoother(dampening);
    	bool enableCoarseCorrection=true;
    	if (enableCoarseCorrection)
      		mgSetup.setupCoarseSuperLUSolver();
    	else
      		mgSetup.setupCoarseNullSolver();
    	mg.levelOperations(levelOp);
    	mg.coarseSolver(mgSetup.coarseSolver());
    	//levelOp.back().maybeRestrictToMaster(newton_rhs);
    	
    	
    	std::function<void(VectorType&)> collect = Dune::ParMG::makeCollect<VectorType>(*mgSetup.comms_.back());
    	std::function<void(VectorType&)> restrictToMaster = [op=levelOp.back()](VectorType& x) { op.maybeRestrictToMaster(x); };
    	//std::cout << "im sweeper vor energyfunctional" << std::endl;
    	

      	auto vz =  A_dune;
	vz *=-1;    	


	auto parallel_energyNorm  = Dune::ParMG::parallelEnergyNorm<VectorType>(vz, restrictToMaster, gridView.grid().comm());




    }*/
              

    template<class EncapsulationTrait>
    size_t
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::size() const
    {
      return this->_size;
    }
  }  // ::pfasst::encap
} // ::pfasst
