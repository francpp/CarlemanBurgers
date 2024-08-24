#include "CarlemanSolver.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include "matrix/kronp.hpp"
#include <iostream>
#include <vector>

namespace sim
{
namespace solvers
{
  // need to pass also the matrice A and dNs
  CarlemanSolver::CarlemanSolver(
    const params::SimulationParameters          &params,
    const discretization::Discretization        &discretization,
    const initial_conditions::InitialConditions &initialConditions)
    : params(params), discretization(discretization),
      initialConditions(initialConditions)
  {}

  void
  CarlemanSolver::solveCarlemanSystem()
  {
    int    nx = params.nx;
    int    nt = params.nt;
    int    N_max = params.N_max;
    double dt = discretization.getTs()[1] - discretization.getTs()[0];
    auto   xs = discretization.getXs();
    auto   ts = discretization.getTs();
  }

} // namespace solvers
} // namespace sim
