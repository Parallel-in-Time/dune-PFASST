#ifndef _ASSEMB_
#define _ASSEMB_

#include <dune/common/function.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include <vector>






template<class LocalView, class MatrixType>
void 

assembleElementA(const LocalView &localView, MatrixType &elementMatrix) {

  using namespace Dune;
  using Element = typename LocalView::Element;

  auto element = localView.element();
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  const auto &localFiniteElement = localView.tree().finiteElement();

  elementMatrix.setSize(localFiniteElement.size(), localFiniteElement.size());
  elementMatrix = 0;




  int order = 2 * (dim * localFiniteElement.localBasis().order() -1 );


  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);


  for (size_t pt = 0; pt < quad.size(); pt++) {   // Loop over all quadrature points


    const auto quadPos = quad[pt].position(); // Position of the current quadrature point in the reference element

    const auto jacobian = geometry.jacobianInverseTransposed(quadPos); // The transposed inverse Jacobian of the map from the reference element to the element

    const auto integrationElement = geometry.integrationElement(quadPos); // The multiplicative factor in the integral transformation formula


    std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos, referenceGradients); // The gradients of the shape functions on the reference element

    std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
    for (size_t i = 0; i < gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);     // Compute the shape function gradients on the real element


    for (size_t i = 0; i < elementMatrix.N(); i++)
      for (size_t j = 0; j < elementMatrix.M(); j++)
        elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
                -= (gradients[i] * gradients[j]) * quad[pt].weight() * integrationElement;
		// += (gradients[i] * gradients[j]) * quad[pt].weight() * integrationElement/(2*PI*PI +1);
  }
}

template<class LocalView, class MatrixType>
void assembleElementM(const LocalView &localView, MatrixType &elementMatrix) {

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
      for (size_t j = 0; j < elementMatrix.M(); j++)
        elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
                += (shapeFunctionValues[i] * shapeFunctionValues[j]) * quad[pt].weight() * integrationElement;

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
                            Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &A,
                            Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &M){
                            //Dune::DenseMatrix<Dune::FieldMatrix<double,1,1>> M_inverse) {

  using namespace Dune;


  auto gridView = basis->gridView();

  MatrixIndexSet occupationPattern;

  getOccupationPattern(basis, occupationPattern);

  occupationPattern.exportIdx(A);
  occupationPattern.exportIdx(M);


  A = 0;
  M = 0;


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


  }
  //M_invers=M;

  //auto isDirichlet = [] (auto x) {return (x[0]<1e-8 or x[0]>0.9999 or x[1]<1e-8 or x[1]>0.9999);}; //ruth_dim

  /*auto isDirichlet = [] (auto x) {return (x[0] < -20.0 + 1e-8 or x[0] > 20.0-1e-8) ;};
  std::vector<char> dirichletNodes;
  interpolate(*basis, dirichletNodes, isDirichlet);  //tu valjda interpoliramo kao na nrpdju

  for (size_t i=0; i<A.N(); i++){
    if (dirichletNodes[i]){
      auto cIt = A[i].begin();
      auto cEndIt = A[i].end();
      for(; cIt!=cEndIt; ++cIt){
        *cIt = (i==cIt.index()) ? 1.0 : 0.0; // 0.0;
      }
    }
  }

  for (size_t i=0; i<M.N(); i++){
    if (dirichletNodes[i]){
      auto cIt = M[i].begin();
      auto cEndIt = M[i].end();
      for(; cIt!=cEndIt; ++cIt){
        *cIt = 0.0;// (i==cIt.index()) ? 1.0 : 0.0;

        //

      }
    }
  }*/

}
#endif
