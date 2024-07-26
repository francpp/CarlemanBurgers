#ifndef MAIN_SIMULATION_HPP
#define MAIN_SIMULATION_HPP

#include "discretization/Discretization.hpp"
#include "error_analysis/ErrorAnalysis.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "matrix/MatrixOperations.hpp"
#include "params/SimulationParameters.hpp"
#include "solvers/CarlemanSolver.hpp"
#include "solvers/EulerSolver.hpp"
#include "solvers/ODE45Solver.hpp"
#include "solvers/PDESolver.hpp"

namespace sim
{
class MainSimulation
{
  params::SimulationParameters          params;
  discretization::Discretization        discretization;
  initial_conditions::InitialConditions initialConditions;
  matrix::MatrixOperations              matrixOperations;

  solvers::EulerSolver    eulerSolver;
  solvers::ODE45Solver    ode45Solver;
  solvers::PDESolver      pdeSolver;
  solvers::CarlemanSolver carlemanSolver;

  error_analysis::ErrorAnalysis errorAnalysis;

public:
  MainSimulation(params::SimulationParameters params);

  void initialize();
  void run();
};
} // namespace sim

#endif // MAIN_SIMULATION_HPP
