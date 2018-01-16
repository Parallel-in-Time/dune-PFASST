 
#ifndef _PFASST__ENCAP__DUNE_VEC_HPP_
#define _PFASST__ENCAP__DUNE_VEC_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;


#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

using Dune::BlockVector;
using Dune::FieldVector;
#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/encap/encapsulation.hpp"

#include <algorithm>

namespace pfasst
{
  namespace encap
  {


  struct dune_encap_tag
          : public encap_data_tag
{};

      /**
      * Spatialized Type Traits for encapsulation of std::vector.
      *
      * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
      * @tparam SpatialPrecision the spatial data precision
      *
      * @ingroup Traits
      */

  template<
          class TimePrecision,
          class SpatialPrecision,
          class gfs, //size_t Dim
	  class M
  >
  struct dune_vec_encap_traits
          //: public encap_traits<TimePrecision, SpatialPrecision, GFS, BlockVector<FieldVector<SpatialPrecision,1>>>; //Dim, BlockVector<FieldVector<SpatialPrecision,1>>>
  {
      using time_t = TimePrecision;
      using spatial_t = SpatialPrecision;

      typedef Dune::YaspGrid<1> Grid;
      typedef Grid::ctype DF;
      typedef Grid::LeafGridView GV;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;

      // Make grid function space
      typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
      typedef Dune::PDELab::istl::VectorBackend<> VBE;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

      typedef double RF; 
      using Z = Dune::PDELab::Backend::Vector<GFS,RF>;	

      using data_t = Z; //BlockVector<FieldVector<spatial_t,1>> ;
      using tag_t = dune_encap_tag;
      using gfs_t = GFS;	
      using mass_t = M;	
      //using dim_t = std::integral_constant<size_t, Dim>;
      //static constexpr size_t  DIM = Dim;
};



/*	//setup the grid
        const int dim=DIM;
	const int degree =1;
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

	typedef double RF; 
  	using Z = Dune::PDELab::Backend::Vector<GFS,RF>;*/




    /**
     * Specialization of Encapsulation for `std::vector`.
     */
    template<
      class EncapsulationTrait
    >
    class Encapsulation<EncapsulationTrait,
                        typename std::enable_if<
                                   std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
                                 >::type>
      :   public std::enable_shared_from_this<Encapsulation<EncapsulationTrait>>
        , public el::Loggable
    {
      public:
        using traits = EncapsulationTrait;
        using factory_t = EncapsulationFactory<traits>;
	//using my_type = decltype(this); Encapsulation<EncapsulationTrait>;

      protected:
        //typename traits::data_t _data;
	shared_ptr<typename traits::data_t> _data;
        const size_t size;

      public:
        explicit Encapsulation(const size_t size = 0);
	explicit Encapsulation(typename EncapsulationTrait::gfs_t gfs);
        Encapsulation(const typename EncapsulationTrait::data_t& data);
        Encapsulation<EncapsulationTrait>& operator=(const typename EncapsulationTrait::data_t& data);

        virtual       typename EncapsulationTrait::data_t& data();
        virtual const typename EncapsulationTrait::data_t& get_data() const;
        virtual size_t get_total_num_dofs() const;
        // assuming square-shaped space
        //virtual std::array<int, EncapsulationTrait::DIM> get_dimwise_num_dofs() const;

        virtual void zero();
        virtual void scaled_add(const typename EncapsulationTrait::time_t& a,
                               const shared_ptr<Encapsulation<EncapsulationTrait>> y);

	virtual void apply_Mass(shared_ptr<typename traits::mass_t> mass, shared_ptr<Encapsulation<EncapsulationTrait>> sol);//, EncapsulationTrait &sol); 

        virtual typename EncapsulationTrait::spatial_t norm0() const;

        template<class CommT>
        bool probe(shared_ptr<CommT> comm, const int src_rank, const int tag);
        template<class CommT>
        void send(shared_ptr<CommT> comm, const int dest_rank, const int tag, const bool blocking);
        template<class CommT>
        void recv(shared_ptr<CommT> comm, const int src_rank, const int tag, const bool blocking);
        template<class CommT>
        void bcast(shared_ptr<CommT> comm, const int root_rank);

        virtual void log(el::base::type::ostream_t& os) const override;
    };

    /**
     * Shortcut for encapsulation of `std::vector` data types.
     */

    template<
      typename time_precision,
      typename spatial_precision,
      typename GFS, //size_t Dim
      typename M	
    >
    using DuneEncapsulation = Encapsulation<dune_vec_encap_traits<time_precision, spatial_precision, GFS, M>>; //Dim>>;


    template<
      class EncapsulationTrait
    >
    class EncapsulationFactory<EncapsulationTrait,
                               typename std::enable_if<
                                          std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
                                        >::type>
      : public std::enable_shared_from_this<EncapsulationFactory<EncapsulationTrait>>
    {
      protected:
        size_t _size;
	typename std::shared_ptr<typename EncapsulationTrait::gfs_t> _gfs;
	//typename EncapsulationTrait::gfs_t *_gfs;

      public:
        explicit EncapsulationFactory(const size_t size = 0);
	explicit EncapsulationFactory(typename EncapsulationTrait::gfs_t gfs);
        EncapsulationFactory(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory(EncapsulationFactory<EncapsulationTrait>&& other);
        virtual ~EncapsulationFactory() = default;
        EncapsulationFactory<EncapsulationTrait>& operator=(const EncapsulationFactory<EncapsulationTrait>& other);
        EncapsulationFactory<EncapsulationTrait>& operator=(EncapsulationFactory<EncapsulationTrait>&& other);

        virtual shared_ptr<Encapsulation<EncapsulationTrait>> create() const;

        virtual void set_size(const size_t& size);
        virtual void set_gfs(typename EncapsulationTrait::gfs_t& gfs);
        virtual size_t size() const;
    };
  }  // ::pfasst::encap
    /*template<typename T>
    static string join(const BlockVector<T>& vec, const string& sep)
    {
#ifndef PFASST_NO_LOGGING
      std::stringstream os;
      os << "[";
      for (size_t i = 0; i < vec.size() - 1; ++i) {
        os << vec[i] << sep;
      }
      os << vec[vec.size()-1];
      os << "]";
      return os.str();
#else
      UNUSED(vec); UNUSED(sep);
    return "";
#endif
    }*/
}  // ::pfasst


#include "pdelab_vec_impl.hpp"

#endif // _PFASST__ENCAP__DUNE_VEC_HPP_
