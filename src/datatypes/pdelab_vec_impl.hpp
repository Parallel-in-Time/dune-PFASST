#include "pdelab_vec.hpp"

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
      : size(size)
    {
        typedef Dune::YaspGrid<1> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<double,1> h = {1};	      
	std::array<int,1> n;
	std::fill(n.begin(), n.end(), size);
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
 	//_data(gfs);
	std::cout << "*************************************   erstelle pdevec mit gfs " << std::endl;
	this->_data = std::make_shared<typename EncapsulationTrait::data_t>(gfs);
		
        //this->zero();
    }


    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(typename EncapsulationTrait::gfs_t gfs)
       : size(100)
    {
	std::cout << "im konstruktor gfs" << std::endl;
      this->_data = std::make_shared<typename EncapsulationTrait::data_t>(gfs);
      //this->zero();
    }

    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(const typename EncapsulationTrait::data_t& data)
   : Encapsulation<EncapsulationTrait>(data.N())
    {
      (this->data()) = data;
    }

    template<class EncapsulationTrait>
    Encapsulation<EncapsulationTrait>&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::operator=(const typename EncapsulationTrait::data_t& data)
    {
      std::cout << "hier im gleich " << std::endl;	
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
	//std::cout << "im data " << std::endl;
	//std::cout << "im data " << _data->N() << std::endl;
      return (*_data);
    }

    template<class EncapsulationTrait>
    const typename EncapsulationTrait::data_t&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_data() const
    {
      return *_data;
    }

    template<class EncapsulationTrait>
    size_t
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_total_num_dofs() const
    {
      return this->get_data().N();
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

      assert(this->get_data().N() == y->data().N());

      std::transform(this->get_data().begin(), this->get_data().end(), y->data().begin(),
                this->data().begin(),
                [a](const typename EncapsulationTrait::spatial_t& xi,
                    const typename EncapsulationTrait::spatial_t& yi) { return xi + a * yi; });
	
	
	/*auto z(*y);  //hier wird nicht wirklich kopiert
	z.data() *= a;
	this->data() += z.data();*/
	
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::spatial_t
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::norm0() const
    {
      return std::abs(*(std::max_element(this->get_data().begin(), this->get_data().end(),
                               [](const typename EncapsulationTrait::spatial_t& a,
                                  const typename EncapsulationTrait::spatial_t& b)
                                 { return std::abs(a) < std::abs(b); })));
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::apply_Mass(shared_ptr<typename traits::mass_t> mass, shared_ptr<Encapsulation<EncapsulationTrait>> sol)//, EncapsulationTrait &sol)
    {
	//std::cout << "matrix vector mult " <<std::endl;
	//std::cout << "matrix vector mult " << Dune::PDELab::Backend::native(this->data())[10][10] <<std::endl;
	//std::cout << "matrix vector mult " << Dune::PDELab::Backend::native((*mass))[1][1][1][1] <<std::endl;
	mass->jacobian_apply((this->data()), (sol->data()));
	//mass.mv(this->data(), sol->data());//sol->data());  
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

    template<class EncapsulationTrait>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::EncapsulationFactory(typename EncapsulationTrait::gfs_t gfs)
      : _gfs(gfs)
    {}


    /*template<class EncapsulationTrait>
    EncapsulationFactory<typename EncapsulationTrait>&  EncapsulationFactory<EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::operator=(const EncapsulationFactory<typename EncapsulationTrait>&& other)
    {	std::cout << "im =" << std::endl;}*/


    /*template<class EncapsulationTrait>
    shared_ptr<Encapsulation<EncapsulationTrait>>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::create() const
    {
	std::cout << "im create" << std::endl;
        typedef Dune::YaspGrid<1> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<double,1> h = {1};	      
	std::array<int,1> n;
	std::cout << "anzahl elemente " << _size << std::endl;
	std::fill(n.begin(), n.end(), _size);
	std::cout << "a " << _size << std::endl;
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(h,n));
	std::cout << "gp erstellt " << _size << std::endl;
        //gridp->refineOptions(false); // keep overlap in cells
        //gridp->globalRefine(1);
        typedef Grid::LeafGridView GV;
	std::cout << "vor gv " << _size << std::endl;
        GV gv=gridp->leafGridView();
	std::cout << "nach gv " << _size << std::endl;
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        FEM fem(gv);

  	// Make grid function space
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	GFS gfs(gv,fem);

	//typedef double RF; 
  	//using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  	//_data(gfs);
		std::cout << "im create vor return " << std::endl;
      return std::make_shared<Encapsulation<EncapsulationTrait>>(gfs);
      //return std::make_shared<EncapsulationTrait>(gfs);	
    }*/


    template<class EncapsulationTrait>
    shared_ptr<Encapsulation<EncapsulationTrait>>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::create() const
    {
	std::cout << "im create" << std::endl;
        
	std::cout << "im create vor return " << std::endl;
      	return std::make_shared<Encapsulation<EncapsulationTrait>>(*(this->_gfs));
      	//return std::make_shared<EncapsulationTrait>(gfs);	
    }

    template<class EncapsulationTrait>
    void
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::set_gfs(typename EncapsulationTrait::gfs_t& gfs)
    {
      this->_size = 100;
      this->_gfs = std::make_shared<typename EncapsulationTrait::gfs_t>(gfs);
    }

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
