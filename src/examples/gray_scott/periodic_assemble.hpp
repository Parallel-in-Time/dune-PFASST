
#include <leathers/push>
#include <leathers/all>

#include <dune/common/function.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/matrixmarket.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include <vector>

#include <leathers/pop>
const size_t NR_OF_COMPONENTS = 2;
//double nu[] = { 2e-5, 1e-5};
double nu[] = { 1e-4, 1e-5};

int inline periodic_pair_right(int node, int size)
{
  
  
  
  return node - size + 1;
}


int inline periodic_pair_up(int node, int size)
{
  
  
  return node - size*(size-1);
}

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

 


  int order = 2 * (dim * localFiniteElement.localBasis().order()-1);


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


  int order = 2 * (dim * localFiniteElement.localBasis().order());

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
  
  int siz = sqrt(basis->size());
   auto isRight = [] (auto x) {return x[0] > 1 - 1e-8;};

    
    auto isUp = [] (auto x) {return x[1] > 1- 1e-8 ;};
   std::vector<int> dirichletLeftNodes, dirichletRightNodes, bottomDirichletNodes, upDirichletNodes ;
    
    interpolate(*basis, dirichletRightNodes, isRight);  //tu valjda interpoliramo kao na nrpdju
  
    interpolate(*basis, upDirichletNodes, isUp);  //tu valjda interpoliramo kao na nrpdju
    
  for (const auto &element : elements(gridView)) { // A loop over all elements of the grid

    localView.bind(element);
    localIndexSet.bind(localView);
    
   
   
  
   
    
    

    for (size_t i = 0; i < localIndexSet.size(); i++) {


      int row = localIndexSet.index(i);
      if(dirichletRightNodes[row])
      {
	 nb.add(row, row);
	 nb.add(row, periodic_pair_right(row, siz));
	 row = periodic_pair_right(row, siz);
      
      }

       if(upDirichletNodes[row] and !dirichletRightNodes[row])
      {
	 nb.add(row, row);
	 nb.add(row, periodic_pair_up(row, siz));
	 row = periodic_pair_up(row, siz);
      
      }
      
      
      for (size_t j = 0; j < localIndexSet.size(); j++) {
	
	
        int col = localIndexSet.index(j);
      
	
	if(dirichletRightNodes[col])
      {  
	 nb.add(col, col);
	 nb.add(col, periodic_pair_right(col, siz));
	 col = periodic_pair_right(col, siz);
	
      
      }

       if(upDirichletNodes[col] and !dirichletRightNodes[col])
      {
	
	
	 nb.add(col, col);
	 nb.add(col, periodic_pair_up(col, siz));
	 col = periodic_pair_up(col, siz);
     
	
      }
	
	
	
	
	
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

 int siz = sqrt(basis->size());
  auto isDirichlet = [] (auto x) {return 0;};
  auto isLeft = [] (auto x) {return x[0] < 1e-8 ;};
  auto isRight = [] (auto x) {return x[0] > 1- 1e-8 ;};

  auto isBottom = [] (auto x) {return x[1] < 1e-8 ;};
  auto isUp = [] (auto x) {return x[1] > 1- 1e-8 ;};

  std::vector<int> dirichletNodes, leftDirichletNodes, rightDirichletNodes, upDirichletNodes, bottomDirichletNodes;

  interpolate(*basis, rightDirichletNodes, isRight);  //tu valjda interpoliramo kao na nrpdju

  interpolate(*basis, upDirichletNodes, isUp);  //tu valjda interpoliramo kao na nrpdju








  for (const auto &element : elements(gridView)) {

    localView.bind(element);
    localIndexSet.bind(localView);
    Matrix <FieldMatrix<double, NR_OF_COMPONENTS,NR_OF_COMPONENTS>> elementMatrix_A;


    assembleElementA(localView, elementMatrix_A);
    Matrix <FieldMatrix<double, NR_OF_COMPONENTS,NR_OF_COMPONENTS>> elementMatrix_M;
    assembleElementM(localView, elementMatrix_M);

    std::vector<int> local_indices;
    local_indices.resize(elementMatrix_A.N());


    for (size_t k = 0; k<elementMatrix_A.N(); k++)
    {
       int row = localIndexSet.index(k);


      if(rightDirichletNodes[row])
	 row=periodic_pair_right(row, siz);



      if(upDirichletNodes[row] and !rightDirichletNodes[row])
	 row=periodic_pair_up(row, siz);

      local_indices[k] = row;

    }



   for (size_t i = 0; i < elementMatrix_A.N(); i++) {

       for (size_t j = 0; j < elementMatrix_A.M(); j++) {

	int row = local_indices[i];
	int col = local_indices[j];
        A[row][col] += elementMatrix_A[i][j];
        M[row][col] += elementMatrix_M[i][j];


      }
    }

  }
  
  /*Dune::FieldMatrix<double,NR_OF_COMPONENTS,NR_OF_COMPONENTS> pom, zero(0.0), neg_pom;
  for(int k=0; k<NR_OF_COMPONENTS; ++k)
  {  pom[k][k] = 1.00*nu[k];
  
    neg_pom[k][k] = -1.00*nu[k];
  }
  for (size_t i=0; i<A.N(); i++){
    if (rightDirichletNodes[i] or upDirichletNodes[i]){
      auto cIt = A[i].begin();
      auto cEndIt = A[i].end();
      
      auto mIt = M[i].begin();
      auto mEndIt = M[i].end();
      for(; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? pom : neg_pom;
      for(; mIt!=mEndIt; ++mIt)
        *mIt = zero;
	
	
   }}*/


  auto isDirichlet2 = [] (auto x) {return (x[0]<1e-8 or x[0]>0.9999 or x[1]<1e-8 or x[1]>0.9999);}; //ruth_dim
  //auto isDirichlet = [] (auto x) {return ( x[0]> 0.99999) or (x[0] < 0.00001) ;};

  std::vector<char> dirichletNodes2;
  interpolate(*basis, dirichletNodes2, isDirichlet2);  //tu valjda interpoliramo kao na nrpdju

  Dune::FieldMatrix<double,NR_OF_COMPONENTS,NR_OF_COMPONENTS> pom, zero(0.0);
  for(int k=0; k<NR_OF_COMPONENTS; ++k)
    pom[k][k] = 1.00;

  for (size_t i=0; i<A.N(); i++){
    if (dirichletNodes2[i]){
      auto cIt = A[i].begin();
      auto cEndIt = A[i].end();
      for(; cIt!=cEndIt; ++cIt){
        *cIt = zero;
      }
    }
  }

  for (size_t i=0; i<M.N(); i++){
    if (dirichletNodes2[i]){
      auto cIt = M[i].begin();
      auto cEndIt = M[i].end();
      for(; cIt!=cEndIt; ++cIt){
        if(i==cIt.index())
          *cIt = pom;
        else
          *cIt = zero;
      }

    }
  }


  /*Dune::FieldMatrix<double,NR_OF_COMPONENTS,NR_OF_COMPONENTS> pom, zero(0.0), neg_pom;
  for(int k=0; k<NR_OF_COMPONENTS; ++k)
  {  pom[k][k] = 1.00*nu[k];

    neg_pom[k][k] = -1.00*nu[k];
  }
  for (size_t i=0; i<A.N(); i++){
    if (rightDirichletNodes[i] or upDirichletNodes[i]){
      auto cIt = A[i].begin();
      auto cEndIt = A[i].end();

      auto mIt = M[i].begin();
      auto mEndIt = M[i].end();
      for(; cIt!=cEndIt; ++cIt)
        *cIt = zero;
      for(; mIt!=mEndIt; ++mIt)
        *mIt = (i==cIt.index()) ? pom : neg_pom;


    }}*/

  
  
 


}