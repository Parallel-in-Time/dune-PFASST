


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


const size_t NR_OF_COMPONENTS = 2;

double nu[] = { 1e-4, 1e-5};

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




  int order = 2 * (dim * localFiniteElement.localBasis().order() );


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

      {

        Dune::FieldMatrix<double,NR_OF_COMPONENTS,NR_OF_COMPONENTS> pom(0.0); // IDENTITY MATRIX
        for (int k = 0; k< NR_OF_COMPONENTS; ++k)
          pom[k][k] =  1.00*nu[k];

        pom*=(gradients[i] * gradients[j]) * quad[pt].weight() * integrationElement;
        elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)]+=pom;

      }
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


  int order = 2 * (dim * localFiniteElement.localBasis().order() );

  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

  for (size_t pt = 0; pt < quad.size(); pt++) {

    const auto quadPos = quad[pt].position();

    const auto integrationElement = geometry.integrationElement(quadPos);


    std::vector<FieldVector<double, 1> > shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    for (size_t i = 0; i < elementMatrix.N(); i++)
      for (size_t j = 0; j < elementMatrix.M(); j++)
      {



        Dune::FieldMatrix<double,NR_OF_COMPONENTS,NR_OF_COMPONENTS> pom(0.0);
        for (int k = 0; k< NR_OF_COMPONENTS; ++k)
          pom[k][k] =  1.00;



        pom*=(shapeFunctionValues[i] * shapeFunctionValues[j]) * quad[pt].weight() * integrationElement;


        elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)] +=pom;

      }
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



template<class Basis, class MatrixType>
void assembleProblem(const Basis &basis,  MatrixType &A, MatrixType &M){
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

    Dune::Matrix <FieldMatrix<double, NR_OF_COMPONENTS,NR_OF_COMPONENTS >> elementMatrix_A;

    //FieldVector<double, nr> mia;
    assembleElementA(localView, elementMatrix_A);
    Dune::Matrix <FieldMatrix<double, NR_OF_COMPONENTS, NR_OF_COMPONENTS>> elementMatrix_M;
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


  /*auto isDirichlet = [] (auto x) {return (x[0]<1e-8 or x[0]>0.9999 or x[1]<1e-8 or x[1]>0.9999);}; //ruth_dim
  //auto isDirichlet = [] (auto x) {return ( x[0]> 0.99999) or (x[0] < 0.00001) ;};

  std::vector<char> dirichletNodes;
  interpolate(*basis, dirichletNodes, isDirichlet);  //tu valjda interpoliramo kao na nrpdju

  Dune::FieldMatrix<double,NR_OF_COMPONENTS,NR_OF_COMPONENTS> pom, zero(0.0);
  for(int k=0; k<NR_OF_COMPONENTS; ++k)
    pom[k][k] = 1.00;

  for (size_t i=0; i<A.N(); i++){
    if (dirichletNodes[i]){
      auto cIt = A[i].begin();
      auto cEndIt = A[i].end();
      for(; cIt!=cEndIt; ++cIt){
        *cIt = zero;
      }
    }
  }

  for (size_t i=0; i<M.N(); i++){
    if (dirichletNodes[i]){
      auto cIt = M[i].begin();
      auto cEndIt = M[i].end();
      for(; cIt!=cEndIt; ++cIt){
        if(i==cIt.index())
          *cIt = pom;
        else
          *cIt = zero;
      }

    }
  }*/


}

