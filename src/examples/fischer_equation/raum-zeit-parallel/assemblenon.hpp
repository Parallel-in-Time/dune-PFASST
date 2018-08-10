


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





template<class LocalView, class MatrixType>
void assembleElementG(const LocalView &localView, VectorType &G, VectorType u) {

  using namespace Dune;
  using Element = typename LocalView::Element;

  auto element = localView.element();
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  const auto &localFiniteElement = localView.tree().finiteElement();


  G.setSize(localFiniteElement.size());
  G = 0;


  int order = 2 * (dim * localFiniteElement.localBasis().order());

  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

  for (size_t pt = 0; pt < quad.size(); pt++) {

    const auto quadPos = quad[pt].position();

    const auto integrationElement = geometry.integrationElement(quadPos);


    std::vector<FieldVector<double, 1> > shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    for (size_t i = 0; i < G.size(); i++)
      for (size_t j = 0; j < G.size(); j++)
        for (size_t k = 0; k < G.size(); k++)
          G[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
                += (shapeFunctionValues[i] * shapeFunctionValues[j] * shapeFunctionValues[k]) * quad[pt].weight() * integrationElement;

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
                            Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &A,
                            Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &M){
                            //Dune::DenseMatrix<Dune::FieldMatrix<double,1,1>> M_inverse) {

  using namespace Dune;


  auto gridView = basis->gridView();

  MatrixIndexSet occupationPattern;

  getOccupationPattern(basis, occupationPattern);

  occupationPattern.exportIdx(A);
  occupationPattern.exportIdx(M);
  occupationPattern.exportIdx(G);

  A = 0;
  M = 0;
  G=0;


  auto localView = basis->localView();
  auto localIndexSet = basis->localIndexSet();



  for (const auto &element : elements(gridView)) {

    localView.bind(element);
    localIndexSet.bind(localView);
    Dune::Matrix <FieldMatrix<double, 1, 1>> elementMatrix_A;
    
    
    assembleElementA(localView, elementMatrix_A);
    Dune::Matrix <FieldMatrix<double, 1, 1>> elementMatrix_M;
    assembleElementM(localView, elementMatrix_M);

    for (size_t i = 0; i < elementMatrix_A.N(); i++) {


      auto row = localIndexSet.index(i);

      for (size_t j = 0; j < elementMatrix_A.M(); j++) {


        auto col = localIndexSet.index(j);


        A[row][col] += elementMatrix_A[i][j];
        M[row][col] += elementMatrix_M[i][j];


      }
    }



    using namespace Dune;
  using Element = typename LocalView::Element;

  auto element = localView.element();
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  const auto &localFiniteElement = localView.tree().finiteElement();


  elementMatrix.setSize(localFiniteElement.size(), localFiniteElement.size());
  elementMatrix = 0;


  int order = 2 * (dim * localFiniteElement.localBasis().order());

  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

  for (size_t pt = 0; pt < quad.size(); pt++) {

    const auto quadPos = quad[pt].position();

    const auto integrationElement = geometry.integrationElement(quadPos);


    std::vector<FieldVector<double, 1> > shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    for (size_t i = 0; i < elementMatrix.N(); i++)
      //for (size_t j = 0; j < elementMatrix.M(); j++)
        elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
                += (shapeFunctionValues[i] * shapeFunctionValues[j]) * quad[pt].weight() * integrationElement;

  }	















  }


}
