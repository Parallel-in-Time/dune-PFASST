// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh>
#include <dune/fufem/assemblers/localassemblers/massassembler.hh>
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh>
#include <dune/fufem/assemblers/istlbackend.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/solvers/iterationsteps/blockgssteps.hh>
#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/common/defaultbitvector.hh>
#include <dune/solvers/common/resize.hh>
#include <dune/solvers/transferoperators/compressedmultigridtransfer.hh>
#include <dune/solvers/iterationsteps/multigridstep.hh>
#include <dune/solvers/solvers/umfpacksolver.hh>
#include <dune/solvers/norms/energynorm.hh>
#include <dune/solvers/norms/twonorm.hh>


#include <dune/tnnmg/iterationsteps/tnnmgstep.hh>
#include <dune/tnnmg/iterationsteps/nonlineargsstep.hh>
#include <dune/tnnmg/functionals/boxconstrainedquadraticfunctional.hh>
#include <dune/tnnmg/functionals/bcqfconstrainedlinearization.hh>
#include <dune/tnnmg/projections/obstacledefectprojection.hh>
#include <dune/tnnmg/localsolvers/scalarobstaclesolver.hh>

// eigene TNNMG Funktionale
#include "examples/FE_Newton/tnnmgfunctional.hh"
#include "examples/FE_Newton/tnnmgfunctionallinearization.hh"
#include "examples/FE_Newton/scalarbisectionsolver.hh"


int main(int argc, char** argv)
{
  // Maybe initialize MPI
  Dune::MPIHelper::instance(argc, argv);


  // setup grid
  constexpr int dim = 2;
  using GridType = Dune::YaspGrid<dim>;
  auto gridptr = Dune::StructuredGridFactory<GridType>::createCubeGrid({0,0},{1,1},{{2,2}});
  gridptr->globalRefine(4); // some refinements

  // setup basis
  using Basis = Dune::Functions::PQ1NodalBasis<GridType::LeafGridView>;
  auto basis = Basis{gridptr->leafGridView()};
  // assemble matrces
  const auto dt = 0.05;
  using MatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1> >;
  MatrixType mass;
  MatrixType stiffness;
  {
    auto sBackend = Dune::Fufem::istlMatrixBackend(stiffness);
    auto mBackend = Dune::Fufem::istlMatrixBackend(mass);

    using Assembler = Dune::Fufem::DuneFunctionsOperatorAssembler<Basis, Basis>;
    auto assembler = Assembler{basis, basis};

    using FiniteElement = std::decay_t<decltype(basis.localView().tree().finiteElement())>;

    auto vintageLaplace = LaplaceAssembler<GridType,FiniteElement, FiniteElement>();
    auto vintageMass = MassAssembler<GridType,FiniteElement, FiniteElement>();

    auto localAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
      vintageLaplace.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
    };
    auto localMassAssembler = [&](const auto& element, auto& localMatrixType, auto&& trialLocalView, auto&& ansatzLocalView){
      vintageMass.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), ansatzLocalView.tree().finiteElement());
      //vintageMass.assemble(element, localMatrixType, trialLocalView.tree().finiteElement(), trialLocalView.tree().finiteElement());
    };

    assembler.assembleBulk(sBackend, localAssembler);
    assembler.assembleBulk(mBackend, localMassAssembler);
  }

  // assemble some (arbitrary) rhs
  using VectorType = Dune::BlockVector<Dune::FieldVector<double, 1> >;
  VectorType rhs(basis.dimension());
  VectorType u(basis.dimension());
  VectorType w(basis.dimension());
  rhs=1.0; // just for testing, TODO: to something sensible
  u=0.0;
  {
    auto tmp=rhs;
    tmp=1;
    mass.mv(tmp, w); // w is just the row sum of mass
  }
  

  // define scalar functions phi, phi', phi''
  const auto _nu = 1.0;
  const auto _n = 1;
  // we do not actually use phi in code
  auto phi = [&](auto&& uu) {return dt*_nu*_nu/(_n+2)*std::pow(uu,_n+2);};

  auto phiprime =
    [&](auto&& uu) {
      return dt*_nu*_nu*std::pow(uu, _n+1);
    };
  auto phi2prime = 
    [&] (auto&& uu) {
      return dt*(_nu*_nu)*(_n+1) * std::pow(uu, _n);
    };


  // setup matrix
  auto matrix_ = mass;
  matrix_ *= (1-dt*_nu*_nu);
  //matrix_.axpy(-dt, stiffness);
  matrix_.axpy(+dt, stiffness); // in contrast to Ruth's code, my stiffness matrix is not negative
 

  // TNNMG
  std::cout << "Start TNNMG Setup" << std::endl;
  using BitVector = Dune::Solvers::DefaultBitVector_t<VectorType>;
  BitVector ignore(u.size());

  //// Transfer setup
  ////
  using TransferOperator = CompressedMultigridTransfer<VectorType>;
  using TransferOperators = std::vector<std::shared_ptr<TransferOperator>>;

  TransferOperators transfer(gridptr->maxLevel());
  for (size_t i = 0; i < transfer.size(); ++i)
  {
    // create transfer operator from level i to i+1
    transfer[i] = std::make_shared<TransferOperator>();
    transfer[i]->setup(*gridptr, i, i+1);
  }

  //// TNNMG without actual constraints
  //using Functional = Dune::TNNMG::EnergyFunctional<MatrixType, VectorType, double,1>;
  using Functional = Dune::TNNMG::EnergyFunctional<MatrixType, VectorType, decltype(phiprime), decltype(phiprime), decltype(phi2prime), double>;

  //auto J = Functional(df, newton_rhs, lower, upper);
  //auto J = Functional(this->M_dune, this->A_dune, dt, _nu, w, rhs);
  auto J = Functional(matrix_, rhs, w, phiprime, phiprime, phi2prime);
  //using LocalFunctional = Dune::TNNMG::EnergyFunctional<MatrixType::block_type, VectorType::block_type, decltype(phiprime), decltype(phiprime), decltype(phi2prime), double>;
  //using LocalFunctional = Dune::TNNMG::EnergyDirectionalRestriction<MatrixType::block_type, VectorType::block_type, decltype(phiprime),double>;

  //auto localSolver = gaussSeidelLocalSolver(Dune::TNNMG::ScalarObstacleSolver());
  //auto localSolver = gaussSeidelLocalSolver(Dune::TNNMG::ScalarBisectionSolver());
  auto localSolver = Dune::TNNMG::ScalarBisectionSolver();

  //auto localSolver = gaussSeidelLocalSolver(TrivialLocalSolver<double,Dune::TNNMG::EnergyFunctional<double,double, DUMMYPHI, decltype(phiprime), decltype(phi2prime), double>, BitVector::value_type>()); // trivial solver does not change anything
  //auto localSolver = TrivialLocalSolver<double,LocalFunctional, BitVector::value_type>(); // trivial solver does not change anything
  //auto localSolver = gaussSeidelLocalSolver(TrivialLocalSolver<double,LocalFunctional, BitVector::value_type>()); // trivial solver does not change anything

  using NonlinearSmoother = Dune::TNNMG::NonlinearGSStep<Functional, decltype(localSolver), BitVector>;
  auto nonlinearSmoother = std::make_shared<NonlinearSmoother>(J, u, localSolver);


  //using Linearization = Dune::TNNMG::BoxConstrainedQuadraticFunctionalConstrainedLinearization<Functional, BitVector>;
  using Linearization = Dune::TNNMG::EnergyFunctionalConstrainedLinearization<Functional, BitVector>; // incorporates what used to be evaluate_f, evaluate_df
  //using DefectProjection = Dune::TNNMG::ObstacleDefectProjection;
  //using LineSearchSolver = TrivialSolver; // always gives 1 as correction damping
  using LineSearchSolver = Dune::TNNMG::ScalarBisectionSolver; // always gives 1 as correction damping
  //using LineSearchSolver = Dune::TNNMG::ScalarObstacleSolver;

  /* Setup linear multigrid */
  using MultiGrid =Dune::Solvers::MultigridStep<MatrixType, VectorType, BitVector>;
  auto mgStep = std::make_shared<MultiGrid>();
  auto gssmoother = Dune::Solvers::BlockGSStepFactory<MatrixType, VectorType, BitVector>::create(Dune::Solvers::BlockGS::LocalSolvers::gs());
  mgStep->setSmoother(&gssmoother);
  mgStep->setTransferOperators(transfer);
  mgStep->setMGType(1,3,3);

  // base solver for multigrid
  auto umfpack = Dune::Solvers::UMFPackSolver<MatrixType, VectorType>{}; // direct solver
  mgStep->basesolver_=&umfpack;

  auto trivialProjection = [](auto&& f, auto& x, auto& c) {};
  //using Step = Dune::TNNMG::TNNMGStep<Functional, BitVector, Linearization, DefectProjection, LineSearchSolver>;
  using Step = Dune::TNNMG::TNNMGStep<Functional, BitVector, Linearization, decltype(trivialProjection), LineSearchSolver>;
  int mu=1;
  //
  // J Functional to minimize, mgStep: linear "solver" for linear correction, mu #multigrid steps per iteration, trivialProjection is identity
  auto step = Step(J, u, nonlinearSmoother, mgStep, mu, trivialProjection, LineSearchSolver());
  step.setIgnore(ignore);

  //using Norm =  TwoNorm<VectorType>;
  using Norm =  EnergyNorm<MatrixType,VectorType>;
  auto norm = Norm(stiffness);
  // max. 20 iterations or two norm of correction less than 1e-10
  using Solver = LoopSolver<VectorType>;
  auto solver = Solver(&step, 20, 1e-10, &norm, Solver::FULL);


  //solver.addCriterion(
  //[&](){
  //return Dune::formatString("   % 12.5e", J(u));
  //},
  //"   energy      ");

  //double initialEnergy = J(u);
  //solver.addCriterion(
  //[&](){
  //static double oldEnergy=initialEnergy;
  //double currentEnergy = J(u);
  //double decrease = currentEnergy - oldEnergy;
  //oldEnergy = currentEnergy;
  //return Dune::formatString("   % 12.5e", decrease);
  //},
  //"   decrease    ");

  solver.addCriterion(
      [&](){
      return Dune::formatString("   % 12.5e", step.lastDampingFactor());
      },
      "   damping     ");


  solver.addCriterion(
      [&](){
      return Dune::formatString("   % 12d", step.linearization().truncated().count());
      },
      "   truncated   ");

  //std::vector<double> correctionNorms;
  //auto tolerance = 1e-8;
  //solver.addCriterion(Dune::Solvers::correctionNormCriterion(step, norm, tolerance, correctionNorms));

  solver.preprocess();
  solver.solve();


  // TNNMG solve
  return 0;
}
