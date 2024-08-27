#ifndef MAIN_SIMULATION_HPP
#define MAIN_SIMULATION_HPP

#include "discretization/Discretization.hpp"
#include "error_analysis/ErrorAnalysis.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include "params/SimulationParameters.hpp"
#include "solvers/CarlemanSolver.hpp"
#include "solvers/EulerSolver.hpp"
#include "solvers/ODE45Solver.hpp"
#include "solvers/PDESolver.hpp"
#include "utils/CarlemanUtils.hpp"
#include "utils/StabilityChecks.hpp"
#include <Eigen/Dense>

namespace sim
{
class MainSimulation
{
  params::SimulationParameters         &params; // Reference to params
  discretization::Discretization        discretization;
  initial_conditions::InitialConditions initialConditions;

  solvers::EulerSolver    eulerSolver;
  solvers::ODE45Solver    ode45Solver;
  solvers::PDESolver      pdeSolver;
  solvers::CarlemanSolver carlemanSolver;

  error_analysis::ErrorAnalysis errorAnalysis;

  Eigen::MatrixXd             F0; // Dense matrix for F0
  Eigen::MatrixXd             F1; // Dense matrix for F1
  Eigen::MatrixXd             F2; // Dense matrix for F2

  void checkStabilityConditions(); // Method to check CFL conditions
  void evaluateCarlemanNumber();   // Method to evaluate Carleman number
  Eigen::MatrixXd convertToDenseEigen(const std::vector<std::vector<double>> &);

public:
  MainSimulation(params::SimulationParameters &params);

  void initialize(); // Method to initialize the simulation
  void run();        // Method to run the simulation
};
} // namespace sim

#endif // MAIN_SIMULATION_HPP
