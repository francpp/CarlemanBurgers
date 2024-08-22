#include "CarlemanSolver.hpp"

namespace sim::solvers
{

CarlemanSolver::CarlemanSolver(
  const params::SimulationParameters               &params,
  const sim::discretization::Discretization        &discretization,
  const sim::initial_conditions::InitialConditions &ic)
{
  // Constructor implementation
}

void
CarlemanSolver::solve()
{
  // Solver implementation
}

void
CarlemanSolver::prepareCarlemanMatrix()
{
  // Prepare Carleman matrix implementation
}

void
CarlemanSolver::solveCarlemanSystem()
{
  // Solve Carleman system implementation
}

} // namespace sim::solvers
