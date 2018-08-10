


#include <dune/common/function.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
//#include <dune/istl/matrixmarket.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include <vector>

const int NR_COMP=1;
typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,NR_COMP,NR_COMP> > MatrixType;
typedef Dune::BlockVector<Dune::FieldVector<double,NR_COMP> > VectorType;



template<class LocalView, class LocalIndexSet>
void assembleElementG(const LocalView &localView, VectorType &G, VectorType &u, LocalIndexSet &localIndexSet) {

  using namespace Dune;
  using Element = typename LocalView::Element;

  auto element = localView.element();
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  const auto &localFiniteElement = localView.tree().finiteElement();

  //LocalIndexSet localIndexSet;

  G.resize(localFiniteElement.size());
  G = 0;


  int order = 2 * (dim * localFiniteElement.localBasis().order());

  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

  for (size_t pt = 0; pt < quad.size(); pt++) {

    const auto quadPos = quad[pt].position();

    const auto integrationElement = geometry.integrationElement(quadPos);


    std::vector<FieldVector<double, 1> > shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    for (size_t i = 0; i < G.size(); i++){
      for (size_t j = 0; j < G.size(); j++){
        for (size_t k = 0; k < G.size(); k++){
        
        
          auto jj = localIndexSet.index(j);
          auto kk = localIndexSet.index(k);
          G[localView.tree().localIndex(i)]
                += (shapeFunctionValues[i] * shapeFunctionValues[j] * u[jj] * shapeFunctionValues[k]) * u[kk] * quad[pt].weight() * integrationElement;
          }}}      

  }
}



template<class Basis>
void getOccupationPattern(const Basis &basis, Dune::MatrixIndexSet &nb) {
  using namespace Dune;
  nb.resize(basis->size(), basis->size());
  auto gridView = basis->gridView();
  auto localView = basis->localView();
  auto localIndexSet = basis->localIndexSet();
  for (const auto &element : elements(gridView)) { // A loop over all elements of the grid

    localView.bind(element);
    localIndexSet.bind(localView);

    for (size_t i = 0; i < localIndexSet.size(); i++) {


      auto row = localIndexSet.index(i);
      for (size_t j = 0; j < localIndexSet.size(); j++) {

        auto col = localIndexSet.index(j);
        nb.add(row, col);
      }
    }
  }

}



template<class Basis>
void assembleProblem(const Basis &basis,
			    VectorType &G,
			    VectorType &u){
                            //Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &A,
                            //Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &M){
                            //Dune::DenseMatrix<Dune::FieldMatrix<double,1,1>> M_inverse) {

  using namespace Dune;


  auto gridView = basis->gridView();

  MatrixIndexSet occupationPattern;

  getOccupationPattern(basis, occupationPattern);

  //occupationPattern.exportIdx(A);
  //occupationPattern.exportIdx(M);
  //occupationPattern.exportIdx(G);

  G.resize(basis->size());

  //A = 0;
  //M = 0;
  G=0;


  auto localView = basis->localView();
  auto localIndexSet = basis->localIndexSet();



  for (const auto &element : elements(gridView)) {

    localView.bind(element);
    localIndexSet.bind(localView);
    //Dune::Matrix <FieldMatrix<double, 1, 1>> elementMatrix_A;
    
    
    //assembleElementA(localView, elementMatrix_A);
    //Dune::Matrix <FieldMatrix<double, 1, 1>> elementMatrix_M;
    //assembleElementM(localView, elementMatrix_M);

    VectorType element_G;	
    assembleElementG(localView, element_G, u, localIndexSet);

    for (size_t i = 0; i < element_G.size(); i++) {



      auto row = localIndexSet.index(i);



      G[row] += element_G[i];	
      //for (size_t j = 0; j < elementMatrix_A.M(); j++) {


        //auto col = localIndexSet.index(j);


        //A[row][col] += elementMatrix_A[i][j];
        //M[row][col] += elementMatrix_M[i][j];


      //}
    }

	





  }


}
