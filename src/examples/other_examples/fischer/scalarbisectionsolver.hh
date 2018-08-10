// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TNNMG_LOCALSOLVERS_SCALARBISECTIONSOLVER_HH
#define DUNE_TNNMG_LOCALSOLVERS_SCALARBISECTIONSOLVER_HH

#include <dune/tnnmg/problem-classes/bisection.hh>


namespace Dune {
namespace TNNMG {



/**
 * \brief A local solver for scalar quadratic obstacle problems
 *
 * \todo Add concept check for the function interface
 */
class ScalarBisectionSolver
{

public:
  ScalarBisectionSolver(double start, double acceptFactor) :
    start_(start),
    acceptFactor_(acceptFactor) {}

  ScalarBisectionSolver() {
    ScalarBisectionSolver(0, 0.9);
  }

  ScalarBisectionSolver(double start) {
    ScalarBisectionSolver(start, 0.9);
  }

  template<class Vector, class Functional, class BitVector>
  void operator()(Vector& x, const Functional& f, const BitVector& ignore) const
  {
    //if (not ignore) // TODO: This still causes errors
    if (true)
    {
      int count=0;
      Bisection bisection(0.0, acceptFactor_);
      x=bisection.minimize(f, start_, 0, count); 
    }
  }
private:
  double start_;
  double acceptFactor_;
};



} // end namespace TNNMG
} // end namespace Dune



#endif // DUNE_TNNMG_LOCALSOLVERS_SCALARBISECTIONSOLVER_HH
