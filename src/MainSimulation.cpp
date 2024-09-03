#include "MainSimulation.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace sim
{

/**
 * @brief Constructs the MainSimulation object and initializes solvers and other
 * components.
 * @param params Reference to the simulation parameters.
 */
MainSimulation::MainSimulation(params::SimulationParameters &params)
  : params(params), // Initialize the member variable params with the reference
    discretization(params), initialConditions(params, discretization),
    eulerSolver(params, discretization, initialConditions),
    ode45Solver(params, discretization, initialConditions),
    pdeSolver(params, discretization, initialConditions),
    carlemanSolver(params, discretization, initialConditions),
    errorAnalysis(eulerSolver, ode45Solver, pdeSolver, carlemanSolver, params)
{}

/**
 * @brief Initializes the simulation, sets up the discretization and initial
 * conditions.
 */
void
MainSimulation::initialize()
{
  params.initialize(); // Initialize params

  // Create the discretization and compute initial conditions
  discretization.createDiscretization();
  checkStabilityConditions();

  initialConditions.computeInitialConditions();
  initialConditions.computeForcingBoundaryConditions(); // Compute the forcing
                                                        // boundary conditions

  // Convert initial conditions to Eigen matrices
  F0 = matrixUtils::convertToDenseEigen(initialConditions.getF0());
  F1 = matrixUtils::convertToDenseEigen(initialConditions.getF1());
  F2 = matrixUtils::convertToDenseEigen(initialConditions.getF2());
}

/**
 * @brief Runs the simulation, solves the systems using different methods, and
 * performs error analysis.
 */
void
MainSimulation::run()
{
  std::cout << "Running simulation..." << std::endl;
  std::cout << params << std::endl;
  std::cout << discretization << std::endl;

  evaluateCarlemanNumber();

  // Solve the Carleman system
  carlemanSolver.solve(F0, F1, F2);

  // Solve the other systems
  eulerSolver.solve(F0, F1, F2);
  ode45Solver.solve(F0, F1, F2);
  pdeSolver.solve(F0, F1, F2);

  // Perform error analysis
  errorAnalysis.computeErrors();

  // Generate plots
  plots::Plotter plotter(params, discretization, initialConditions,
                         carlemanSolver, eulerSolver, pdeSolver, ode45Solver,
                         errorAnalysis);
  plotter.initialize();
  plotter.plotSolution();
  plotter.plotErrors();
  plotter.plotErrorConvergence();
}

/**
 * @brief Checks the stability conditions, such as CFL conditions, for the
 * simulation.
 * @throws std::runtime_error if the CFL conditions are not met.
 */
void
MainSimulation::checkStabilityConditions()
{
  try
    {
      sim::stability::checkCFLConditions(
        params.U0,
        discretization.getTs()[1] - discretization.getTs()[0], // dt
        discretization.getXs()[1] - discretization.getXs()[0], // dx
        params.nu,
        (params.T) / (params.nt * 10 - 1), // dt_ode
        (params.L0) / (params.nx_pde - 1), // dx_pde
        (params.T) / (params.nt_pde - 1)   // dt_pde
      );
    }
  catch(const std::runtime_error &e)
    {
      std::cerr << "CFL Condition Error: " << e.what() << std::endl;
      throw; // Re-throw the exception to terminate the simulation
    }
}

/**
 * @brief Evaluates the Carleman convergence number.
 * @throws std::runtime_error if there is an error in the calculation.
 */
void
MainSimulation::evaluateCarlemanNumber()
{
  try
    {
      double R = sim::utils::calculateCarlemanConvergenceNumber(
        F0, F1, F2, initialConditions.getU0s(),
        discretization.getTs()[1] - discretization.getTs()[0], // dt
        params.nt, params.nx, params.N_max);
      std::cout << "Carleman convergence number R: " << R << std::endl;
    }
  catch(const std::runtime_error &e)
    {
      std::cerr << "Carleman Convergence Number Error: " << e.what()
                << std::endl;
      throw; // Re-throw the exception to terminate the simulation
    }
}

} // namespace sim
