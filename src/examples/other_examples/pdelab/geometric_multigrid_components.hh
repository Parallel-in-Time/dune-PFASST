#ifndef DUNE_PARSOLVE_GEOMETRICMULTIGRIDCOMPONENTS_HH
#define DUNE_PARSOLVE_GEOMETRICMULTIGRIDCOMPONENTS_HH

#include<utility>
#include<vector>
#include<set>
#include<map>

#include <dune/pdelab/boilerplate/pdelab.hh>


/** Hierarchy of continuous Lagrange finite element spaces

    Constructs a grid function space for all levels based on level grid views.
    Setting the current level allows the use like a grid function
    space constructed on the leaf view without giving further parameters.

    \param T Grid
    \param N number type
*/
/*template<typename T, typename N, unsigned int degree, typename BCType,
         Dune::GeometryType::BasicType gt,
         Dune::PDELab::MeshType mt, Dune::SolverCategory::Category st = Dune::SolverCategory::sequential,
         typename VBET=Dune::PDELab::istl::VectorBackend<> >
class CGSpaceHierarchy {
public:
  // export types
  typedef T Grid;
  typedef typename T::LevelGridView GV;
  typedef typename T::ctype ctype;
  static const int dim = T::dimension;
  static const int dimworld = T::dimensionworld;

  typedef Dune::PDELab::CGFEMBase<GV,ctype,N,degree,dim,gt> FEMB;
  typedef Dune::PDELab::CGCONBase<Grid,degree,gt,mt,st,BCType> CONB;

  typedef typename FEMB::FEM FEM;
  typedef typename CONB::CON CON;

  typedef VBET VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

  typedef N NT;
  using DOF = Dune::PDELab::Backend::Vector<GFS,N>;
  typedef Dune::PDELab::DiscreteGridFunction<GFS,DOF> DGF;
  typedef typename GFS::template ConstraintsContainer<N>::Type CC;
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> VTKF;

  // constructor making the grid function space an all that is needed
  CGSpaceHierarchy (Grid& grid_, const BCType& bctype)
    : grid(grid_), gvp(grid_.maxLevel()+1), fembp(grid_.maxLevel()+1), conbp(grid_.maxLevel()+1),
      gfsp(grid_.maxLevel()+1), ccp(grid_.maxLevel()+1)
  {
    for (int j=0; j<=grid.maxLevel(); j++)
      {
        gvp[j] = Dune::shared_ptr<GV>(new GV(grid.levelGridView(j)));
        fembp[j] = Dune::shared_ptr<FEMB>(new FEMB(*gvp[j]));
        conbp[j] = Dune::shared_ptr<CONB>(new CONB(grid,bctype));
        gfsp[j] = Dune::shared_ptr<GFS>(new GFS(*gvp[j],fembp[j]->getFEM(),conbp[j]->getCON()));
        conbp[j]->postGFSHook(*gfsp[j]);
        ccp[j] = Dune::shared_ptr<CC>(new CC());
      }
  }

  // get max level
  int maxLevel () const
  {
    return grid.maxLevel();
  }

  // return finite element map reference
  FEM& getFEM(int j)
  {
    return fembp[j].getFEM();
  }

  // return finite element map reference const version
  const FEM& getFEM(int j) const
  {
    return fembp[j].getFEM();
  }

  // return gfs reference
  GFS& getGFS (int j)
  {
    return *gfsp[j];
  }

  // return gfs reference const version
  const GFS& getGFS (int j) const
  {
    return *gfsp[j];
  }

  // return gfs reference
  CC& getCC (int j)
  {
    return *ccp[j];
  }

  // return gfs reference const version
  const CC& getCC (int j) const
  {
    return *ccp[j];
  }

  // assemble constraints for level j
  void assembleConstraints (int j, const BCType& bctype)
  {
    ccp[j]->clear();
    Dune::PDELab::constraints(bctype,*gfsp[j],*ccp[j]);
  }

  // assemble all constraints for all levels
  void assembleAllConstraints (const BCType& bctype)
  {
    for (int j=0; j<=grid.maxLevel(); j++)
      {
        ccp[j]->clear();
        Dune::PDELab::constraints(bctype,*gfsp[j],*ccp[j]);
      }
  }

  // reset constraints
  void clearConstraints (int j)
  {
    ccp[j]->clear();
  }

  // reset all constraints
  void clearAllConstraints ()
  {
    for (int j=0; j<=grid.maxLevel(); j++)
      ccp[j]->clear();
  }

  void setConstrainedDOFS (int j, DOF& x, NT nt) const
  {
    Dune::PDELab::set_constrained_dofs(*ccp[j],nt,x);
    conbp[j]->make_consistent(*gfsp[j],x);
  }

  void setNonConstrainedDOFS (int j, DOF& x, NT nt) const
  {
    Dune::PDELab::set_nonconstrained_dofs(*ccp[j],nt,x);
    conbp[j]->make_consistent(*gfsp[j],x);
  }

  void copyConstrainedDOFS (int j, const DOF& xin, DOF& xout) const
  {
    Dune::PDELab::copy_constrained_dofs(*ccp[j],xin,xout);
    conbp[j]->make_consistent(*gfsp[j],xout);
  }

  void copyNonConstrainedDOFS (int j, const DOF& xin, DOF& xout) const
  {
    Dune::PDELab::copy_nonconstrained_dofs(*ccp[j],xin,xout);
    conbp[j]->make_consistent(*gfsp[j],xout);
  }

private:
  Grid& grid;
  // vectors with shared pointers gridview, fem, constraints,gfs etc.
  std::vector<Dune::shared_ptr<GV> > gvp;
  std::vector<Dune::shared_ptr<FEMB> > fembp;
  std::vector<Dune::shared_ptr<CONB> > conbp;
  std::vector<Dune::shared_ptr<GFS> > gfsp;
  std::vector<Dune::shared_ptr<CC> > ccp;
  // make copy constructor private to disallow copying
  CGSpaceHierarchy (const CGSpaceHierarchy& other) {}
};*/



/** multigrid transfer operator assembled as a matrix

    Uses a grid function space hierarchy to construct prolongation matrix
    between two CONSECUTIVE levels of the hierarchy. Works for arbitrary
    polynomial degree and spatial dimension.
*/
template<typename GFS>
class ProlongationOperator : public Dune::BCRSMatrix< Dune::FieldMatrix< double,1,1> >
{
  // export types
  typedef double E;
  typedef Dune::FieldMatrix<E,1,1> M;
  //typedef typename GFSH::GFS GFS;
  typedef Dune::YaspGrid<1> Grid;
  typedef typename Grid::LevelGridView GV;
  typedef typename GV::Grid::ctype ctype;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::BCRSMatrix<M>::size_type size_type;
  typedef std::set<size_type> IndexSet;
  typedef std::vector<IndexSet> Graph;

  class CoordinateEvaluation
  {
  public:
    // store the coordinate to evaluate
    CoordinateEvaluation (int i_) : i(i_) {}

    // eval coordinate i
    template<typename DT, typename RT>
    inline void evaluate (const DT& x, RT& y) const
    {
      y = x[i];
      return;
    }

  private:
    int i;
  };

public:

  typedef E ElementType;
  typedef Dune::BCRSMatrix<M> BaseT;

  ProlongationOperator (const GFS& gfs_1, const GFS& gfs_0)
    : BaseT(gfs_1.globalSize(),gfs_0.globalSize(),
            Dune::BCRSMatrix<M>::random),
      gfsc(gfs_1), gfsf(gfs_0)
  {
    // check level
    //if (level<=0 || level>gfsh.maxLevel())
    //  DUNE_THROW(Dune::Exception, "level out of range");

    // gfs and matrix graph
    // Dune::shared_ptr<const GFS> pgfsf(gfsh_.getGFS(level));
    // Dune::shared_ptr<const GFS> pgfsc(gfsh_.getGFS(level-1));
    //const GFS& gfsf = gfsh_.getGFS(level);  // fine gfs
    //const GFS& gfsc = gfsh_.getGFS(level-1);// coarse gfs
    Graph graph(gfsf.globalSize());                        // matrix graph

    // make a vector on the fine grid containing the global indices
    using X = Dune::PDELab::Backend::Vector<GFS,E>;
    using Dune::PDELab::Backend::native;
    X xf(gfsf,0.0);
    for (typename X::size_type i=0; i<xf.N(); i++) native(xf)[i] = i;
    X xc(gfsc,0.0);
    for (typename X::size_type i=0; i<xc.N(); i++) native(xc)[i] = i;

    // make local function spaces
    typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
    LFS lfsf(gfsf);
    LFS lfsc(gfsc);
    typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
    LFSCache lfsf_cache(lfsf);
    LFSCache lfsc_cache(lfsc);
    typedef typename X::template LocalView<LFSCache> LView;
    LView lviewf(xf);
    LView lviewc(xc);
    std::cout << "vor der ersten schlefe" << std::endl;	
    // loop over fine grid to get matrix graph
    for (const auto& cell : elements(gfsf.gridView()))
      {
        // bind local function spaces to element in fine and coarse grid
        lfsf.bind(cell);
        lfsf_cache.update();
        lfsc.bind(cell.father());
        lfsc_cache.update();

        // get global indices from the helper vector
        std::vector<E> indexf(lfsf.size());
        lviewf.bind(lfsf_cache);
        lviewf.read(indexf);
        lviewf.unbind();
        std::vector<E> indexc(lfsc.size());
        lviewc.bind(lfsc_cache);
        lviewc.read(indexc);
        lviewc.unbind();
        std::cout << "determine local coordinates of vertices in fine element" << std::endl;	
        // determine local coordinates of vertices in fine element
        const int dim = GV::dimension;
        std::vector<Dune::FieldVector<ctype,dim> > local_position(lfsf.size());
        for (int k=0; k<dim; k++)
          {
            CoordinateEvaluation f(k);
            std::vector<ctype> c(lfsf.size());
            lfsf.finiteElement().localInterpolation().interpolate(f,c);
            for (typename LFS::Traits::SizeType i=0; i<lfsf.size(); i++) local_position[i][k] = c[i];
          }

        // determine local coordinates of vertices in father element
        std::vector<Dune::FieldVector<ctype,dim> > local_position_in_father(lfsf.size());
        for (typename LFS::Traits::SizeType i=0; i<lfsf.size(); i++)
          local_position_in_father[i] = cell.geometryInFather().global(local_position[i]);

        // interpolation weights are values of coarse grid basis functions at fine grid points
        for (typename LFS::Traits::SizeType i=0; i<lfsf.size(); i++)
          {
            typedef typename LFS::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::RangeType RangeType;
            std::vector<RangeType> phi(lfsc.size());
            lfsc.finiteElement().localBasis().evaluateFunction(local_position_in_father[i],phi);
            for (typename LFS::Traits::SizeType j=0; j<lfsc.size(); j++)
              {
                if (phi[j]>1E-6)
                  graph[static_cast<size_type>(indexf[i])].insert(static_cast<size_type>(indexc[j]));
              }
          }
      }
    std::cout << "    // now set up the sparse matrix pattern" << std::endl;	
	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
    // now set up the sparse matrix pattern
    for (typename Graph::size_type i=0; i<graph.size(); ++i){
      this->setrowsize(i,graph[i].size()); std::cout << " wichtig  " <<  i  << " " <<  graph[i].size() << std::endl; }
    this->endrowsizes();     std::cout << "   endsetrow "  <<  graph.size() << std::endl;	
    for (typename Graph::size_type i=0; i<graph.size(); ++i)
      {
	std::cout << " +++++++++++++++++++++++++++++++++++++++++ start " <<   std::endl;
	graph[0].end();
	typename IndexSet::iterator it=graph[1].begin(); ++it;

	std::cout << " +++++++++++++++++++++++++++++++++++++++++ end " <<   std::endl;
	//if(*graph[0].begin() == *graph[0].end())  	  std::cout << " +++++++++++++++++++++++++++++++++++++++++ sollte ok sein " <<   std::endl;
	int j=0;
        for (typename IndexSet::iterator it=graph[i].begin(); it!=graph[i].end();  ++it){
	  std::cout << my_rank << " vor  enddindex " <<  (double) i << std::endl;
          this->addindex(i,*it);   std::cout << "   enddindex " << std::endl;
	 std::cout << "   nach innerer schleife "  << std::endl;
	  if(i%2 ==0)break; if(j==1)break;  j++; 
	}
	if(i==1000) break;	
      } std::cout << my_rank << "   nach schleife "  << std::endl; //std::exit(0);
    this->endindices();
    std::cout << my_rank << "  beide noch am start mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm  loop over grid for the second time and insert values " << std::endl;	
    // loop over grid for the second time and insert values
    int ruth = 0;	
    for (const auto& cell : elements(gfsf.gridView())) //problem loop
      {
	std::cout << my_rank << " ------------------------------------------------------------------------------------------------  ruth " << ruth << std::endl; if(ruth==1000) break;  ruth++;
        // bind local function spaces to element in fine and coarse grid
        lfsf.bind(cell);
        lfsf_cache.update();
        lfsc.bind(cell.father());
        lfsc_cache.update();

        // get global indices from the helper vector
        std::vector<E> indexf(lfsf.size());
        lviewf.bind(lfsf_cache);
        lviewf.read(indexf);
        lviewf.unbind();
        std::vector<E> indexc(lfsc.size());
        lviewc.bind(lfsc_cache);
        lviewc.read(indexc);
        lviewc.unbind();

        // determine local coordinates of vertices in fine element
        const int dim = GV::dimension;
        std::vector<Dune::FieldVector<ctype,dim> > local_position(lfsf.size());
        for (int k=0; k<dim; k++)
          {
            CoordinateEvaluation f(k);
            std::vector<ctype> c(lfsf.size());
            lfsf.finiteElement().localInterpolation().interpolate(f,c);
            for (typename LFS::Traits::SizeType i=0; i<lfsf.size(); i++)
              local_position[i][k] = c[i];
          }
        std::cout << "   determine local coordinates of vertices in father element " << std::endl;	
        // determine local coordinates of vertices in father element
        std::vector<Dune::FieldVector<ctype,dim> > local_position_in_father(lfsf.size());
        for (typename LFS::Traits::SizeType i=0; i<lfsf.size(); i++)
          local_position_in_father[i] = cell.geometryInFather().global(local_position[i]);
        std::cout << my_rank << "  interpolate " << std::endl;	
        // interpolation weights are values of coarse grid basis functions at fine grid points
        for (typename LFS::Traits::SizeType i=0; i<lfsf.size(); i++)
          {
            typedef typename LFS::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::RangeType RangeType;
            std::vector<RangeType> phi(lfsc.size());
            lfsc.finiteElement().localBasis().evaluateFunction(local_position_in_father[i],phi);
            for (typename LFS::Traits::SizeType j=0; j<lfsc.size(); j++)
              if (phi[j]>1E-6){
		std::cout << my_rank << " vor dieser zuweisung " << i << " " << j << std::endl;
                (*this)[static_cast<size_type>(indexf[i])][static_cast<size_type>(indexc[j])] = (double) phi[j]; std::cout << my_rank << " nach dieser zuweisung " << i << " " << j << std::endl;}  //das hier macht die probleme!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          }        std::cout << my_rank << "  interpolate " << std::endl;	
      }
	std::cout << my_rank << " ende konstruktor " <<  std::endl;
  }//ende konstruktor prolongation operator

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ProlongationOperator& operator= (const E& x)
  {
    BaseT::operator=(x);
    return *this;
  }

  // for debugging and AMG access
  BaseT& base ()
  {
    return *this;
  }

  const BaseT& base () const
  {
    return *this;
  }

private:
  const GFS& gfsc;
  const GFS& gfsf;
  //int level;
};


/** Manage a hierarchy of prolongation operators for a given grid function space hierarchy
 */
/*template<typename GFSH>
class ProlongationOperatorHierarchy
{
public:

  // prolongation operator type
  typedef ProlongationOperator<GFSH> PO;

  ProlongationOperatorHierarchy (const GFSH& gfsh_)
    : gfsh(gfsh_), pmp(gfsh_.maxLevel()+1)
  {
    for (int j=1; j<=gfsh.maxLevel(); j++)
      pmp[j] = Dune::shared_ptr<PO>(new PO(gfsh,j));
  }

  // get a prolongation operator
  const PO& getPO (int level) const
  {
    if (level>=1 && level<pmp.size())
      return *pmp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a prolongation operator
  const PO& operator[] (int level) const
  {
    if (level>=1 && level<pmp.size())
      return *pmp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a prolongation operator
  PO& getPO (int level)
  {
    if (level>=1 && level<pmp.size())
      return *pmp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a prolongation operator
  PO& operator[] (int level)
  {
    if (level>=1 && level<pmp.size())
      return *pmp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  int maxLevel () const
  {
    return gfsh.maxLevel();
  }

private:
  const GFSH& gfsh;
  std::vector<Dune::shared_ptr<PO> > pmp;
  // make copy constructor private to disallow copying
  ProlongationOperatorHierarchy (const ProlongationOperatorHierarchy& other) {}
};*/



/**
   Manage a hierarchy of degree of freedom vectors for a given grid function space hierarchy
*/
/*template<typename GFSH>
class VectorHierarchy
{
  // constants and types
  typedef typename GFSH::GFS GFS;

public:
  // constraints container type
  typedef typename GFSH::DOF Vector;

  VectorHierarchy (const GFSH& gfsh_) : gfsh(gfsh_), vectorp(gfsh_.maxLevel()+1)
  {
    for (int j=0; j<=gfsh.maxLevel(); j++)
      vectorp[j] = Dune::shared_ptr<Vector>(new Vector(gfsh.getGFS(j)));
  }

  // get a const vector
  const Vector& getVector (int level) const
  {
    if (level>=0 && level<vectorp.size())
      return *vectorp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a vector
  Vector& getVector (int level)
  {
    if (level>=0 && level<vectorp.size())
      return *vectorp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a const vector
  const Vector& operator[] (int level) const
  {
    if (level>=0 && level<vectorp.size())
      return *vectorp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a vector
  Vector& operator[] (int level)
  {
    if (level>=0 && level<vectorp.size())
      return *vectorp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  int maxLevel () const
  {
    return gfsh.maxLevel();
  }

private:
  const GFSH& gfsh;
  std::vector<Dune::shared_ptr<Vector> > vectorp;
  // make copy constructor private to disallow copying
  VectorHierarchy (const VectorHierarchy& other) {}
};*/


/**
   Manage a hierarchy of stiffness matrices for a given grid function space hierarchy
*/
/*template<typename GOSH>
class MatrixHierarchy
{
public:
  typedef typename GOSH::GO GO;
  typedef typename GOSH::MAT Matrix;


  MatrixHierarchy (const GOSH& gosh_)
    : gosh(gosh_), matrixp(gosh.maxLevel()+1)
  {
    for (int l=0; l<=gosh.maxLevel(); l++)
      matrixp[l] = Dune::shared_ptr<Matrix>(new Matrix(gosh.getGO(l)));
  }

  // get a const matrix
  const Matrix& getMatrix (int level) const
  {
    if (level>=0 && level<matrixp.size())
      return *matrixp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a matrix
  Matrix& getMatrix (int level)
  {
    if (level>=0 && level<matrixp.size())
      return *matrixp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a const matrix
  const Matrix& operator[] (int level) const
  {
    if (level>=0 && level<matrixp.size())
      return *matrixp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // get a matrix
  Matrix& operator[] (int level)
  {
    if (level>=0 && level<matrixp.size())
      return *matrixp[level];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  int maxLevel () const
  {
    return gosh.maxLevel();
  }

  bool check () const
  {
    for (int l=0; l<matrixp.size(); l++)
      {
        for (int i=0; i<matrixp[l]->N(); i++)
          if ((*matrixp[l])[i][i]>10.0)
            {
              std::cout << "found matrix element larger than 10, level=" << l << " index=" << i << " " << (*matrixp[l])[i][i] << std::endl;
              return true;
            }
      }
    return false;
  }

private:
  const GOSH& gosh;
  std::vector<Dune::shared_ptr<Matrix> > matrixp;

  // make copy constructor private to disallow copying
  MatrixHierarchy (const MatrixHierarchy& other) {}
};*/


/**
   Manage a hierarchy of assemblers
*/
/*template<typename FS, typename LOP, Dune::SolverCategory::Category st = Dune::SolverCategory::sequential>
class GalerkinGlobalAssemblerHierarchy
{
public:
  // export types
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  typedef Dune::PDELab::GridOperator<typename FS::GFS,typename FS::GFS,LOP,MBE,
                                     typename FS::NT,typename FS::NT,typename FS::NT,
                                     typename FS::CC,typename FS::CC> GO;
  typedef typename GO::Jacobian MAT;
  typedef VectorHierarchy<FS> VH;
  typedef MatrixHierarchy<GalerkinGlobalAssemblerHierarchy<FS,LOP,st> > MH;

  GalerkinGlobalAssemblerHierarchy (const FS& fs, LOP& lop, const std::size_t nonzeros) : gop(fs.maxLevel()+1)
  {
    for (int j=0; j<=fs.maxLevel(); j++)
      gop[j] = Dune::shared_ptr<GO>(new GO(fs.getGFS(j),fs.getCC(j),fs.getGFS(j),fs.getCC(j),lop,MBE(nonzeros)));
  }

  int maxLevel () const
  {
    return gop.size()-1;
  }

  // return grid reference
  GO& getGO (int j)
  {
    if (j>=0 && j<gop.size())
      return *gop[j];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // return grid reference const version
  const GO& getGO (int j) const
  {
    if (j>=0 && j<gop.size())
      return *gop[j];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // evaluate jacobian on all levels
  void jacobian (const VH& xh, MH& mh)
  {
    for (int j=0; j<=maxLevel(); j++)
      {
        mh[j] = 0.0;
        gop[j]->jacobian(xh[j],mh[j]);
      }
  }

private:
  std::vector<Dune::shared_ptr<GO> > gop;
};

/** Manage a hierarchy of assemblers (nonoverlapping variant)
 */
/*template<typename FS, typename LOP>
class GalerkinGlobalAssemblerHierarchy<FS,LOP,Dune::SolverCategory::nonoverlapping>
{
public:
  // export types
  typedef typename FS::VBE::MatrixBackend MBE;
  typedef Dune::PDELab::GridOperator<typename FS::GFS,typename FS::GFS,LOP,MBE,
                                     typename FS::NT,typename FS::NT,typename FS::NT,
                                     typename FS::CC,typename FS::CC> GO; //ruth 

  typedef typename GO::Jacobian MAT;
  typedef VectorHierarchy<FS> VH;
  typedef MatrixHierarchy<GalerkinGlobalAssemblerHierarchy<FS,LOP,Dune::SolverCategory::nonoverlapping> > MH;

  GalerkinGlobalAssemblerHierarchy (const FS& fs, LOP& lop) : gop(fs.maxLevel()+1)
  {
    for (int j=0; j<=fs.maxLevel(); j++)
      gop[j] = Dune::shared_ptr<GO>(new GO(fs.getGFS(j),fs.getCC(j),fs.getGFS(j),fs.getCC(j),lop));
  }

  int maxLevel () const
  {
    return gop.size()-1;
  }

  // return grid reference
  GO& getGO (int j)
  {
    if (j>=0 && j<gop.size())
      return *gop[j];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // return grid operator reference const version
  const GO& getGO (int j) const
  {
    if (j>=0 && j<gop.size())
      return *gop[j];
    else
      DUNE_THROW(Dune::Exception, "level out of range");
  }

  // evaluate jacobian on all levels
  void jacobian (const VH& xh, MH& mh)
  {
    for (int j=0; j<=maxLevel(); j++)
      {
        mh[j] = 0.0;
        gop[j]->jacobian(xh[j],mh[j]);
      }
  }

private:
  std::vector<Dune::shared_ptr<GO> > gop;
};*/

#endif

