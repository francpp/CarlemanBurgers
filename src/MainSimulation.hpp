#ifndef MAIN_SIMULATION_HPP
#define MAIN_SIMULATION_HPP

#include "discretization/Discretization.hpp"
#include "error_analysis/ErrorAnalysis.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include "params/SimulationParameters.hpp"
#include "plots/plots.hpp"
#include "solvers/CarlemanSolver.hpp"
#include "solvers/EulerSolver.hpp"
#include "solvers/ODE45Solver.hpp"
#include "solvers/PDESolver.hpp"
#include "utils/CarlemanUtils.hpp"
#include "utils/MatrixFormat.hpp"
#include "utils/StabilityChecks.hpp"

#include <Eigen/Dense>

namespace sim
{
/**
 * @class MainSimulation
 * @brief This class manages the entire simulation process, including
 * initialization, running the simulation, and error analysis.
 */
class MainSimulation
{
  params::SimulationParameters &params; ///< Reference to simulation parameters.
  discretization::Discretization
    discretization; ///< Handles the discretization of the domain.
  initial_conditions::InitialConditions
    initialConditions; ///< Manages the initial and boundary conditions.

  solvers::EulerSolver eulerSolver; ///< Solver using Euler's method.
  solvers::ODE45Solver ode45Solver; ///< Solver using MATLAB's ODE45 method.
  solvers::PDESolver
    pdeSolver; ///< Solver for the PDE using central difference.
  solvers::CarlemanSolver
    carlemanSolver; ///< Solver using the Carleman linearization method.

  error_analysis::ErrorAnalysis
    errorAnalysis; ///< Performs error analysis on the results.

  Eigen::MatrixXd F0; ///< Dense matrix for initial conditions.
  Eigen::MatrixXd F1; ///< Dense matrix for first-order conditions.
  Eigen::MatrixXd F2; ///< Dense matrix for second-order conditions.

  /**
   * @brief Checks the stability conditions, such as CFL conditions, for the
   * simulation.
   * @throws std::runtime_error if the CFL conditions are not met.
   */
  void checkStabilityConditions();

  /**
   * @brief Evaluates the Carleman convergence number.
   * @throws std::runtime_error if there is an error in the calculation.
   */
  void evaluateCarlemanNumber();

public:
  /**
   * @brief Constructor for MainSimulation.
   * @param params Reference to the simulation parameters.
   */
  MainSimulation(params::SimulationParameters &params);

  /**
   * @brief Initializes the simulation by setting up parameters, discretization,
   * and initial conditions.
   */
  void initialize();

  /**
   * @brief Runs the simulation, solves the systems, and performs error
   * analysis.
   */
  void run();
};

} // namespace sim

#endif // MAIN_SIMULATION_HPP
