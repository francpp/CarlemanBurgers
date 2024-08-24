#include "MainSimulation.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include <iostream>

namespace sim
{

MainSimulation::MainSimulation(params::SimulationParameters &params)
  : params(params), // Initialize the member variable params with the reference
    discretization(params), initialConditions(params, discretization),
    eulerSolver(params, discretization, initialConditions),
    ode45Solver(params, discretization, initialConditions),
    pdeSolver(params, discretization, initialConditions),
    carlemanSolver(params, discretization, initialConditions),
    errorAnalysis(eulerSolver, ode45Solver, pdeSolver, carlemanSolver)
{
  // Constructor implementation
}

void
MainSimulation::initialize()
{
  params
    .initialize(); // Initialize params (which will now affect all components)

  // Create the discretization and compute initial conditions
  discretization.createDiscretization();
  initialConditions.computeInitialConditions();
  initialConditions.computeForcingBoundaryConditions(); // Compute the forcing
                                                        // boundary conditions

  // Initialize F0, F1, F2 using the initial conditions
  /*F0 = Eigen::SparseMatrix<double>(params.nt, params.nx);
  F1 = Eigen::MatrixXd::Zero(params.nx, params.nx);
  F2 = Eigen::MatrixXd::Zero(params.nx, params.nx * params.nx);

  // Fill F0, F1, F2 with data from initialConditions
  for (int i = 0; i < params.nt; ++i) {
    for (int j = 0; j < params.nx; ++j) {
      F0.insert(i, j) = initialConditions.getF0()[i][j];
    }
  }

  for (int i = 0; i < params.nx; ++i) {
    for (int j = 0; j < params.nx; ++j) {
      F1(i, j) = initialConditions.getF1()[i][j];
    }
    for (int j = 0; j < params.nx * params.nx; ++j) {
      F2(i, j) = initialConditions.getF2()[i][j];
    }
  }*/
}

Eigen::SparseMatrix<double>
MainSimulation::prepareCarlemanMatrix()
{
  std::vector<int> dNs = matrix::calculateBlockSizes(params.N_max, params.nx);
  return matrix::assembleCarlemanMatrix(dNs, params.N_max, params.nx,
                                        params.ode_deg, F0, F1, F2);
}

void
MainSimulation::run()
{
  std::cout << "Running simulation..." << std::endl;
  std::cout << params << std::endl;
  std::cout << discretization << std::endl;
  std::cout << initialConditions << std::endl;

  checkStabilityConditions();

  // Prepare the Carleman matrix
  Eigen::SparseMatrix<double> carlemanMatrix = prepareCarlemanMatrix();

  evaluateCarlemanNumber();

  // Solve the Carleman system
  carlemanSolver.solveCarlemanSystem();

  // Proceed with the rest of your simulation process...
}

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
        (params.T) / (params.nt * 10 - 1),     // dt_ode
        (params.L0 / 2) / (params.nx_pde - 1), // dx_pde
        (params.T) / (params.nt_pde - 1)       // dt_pde
      );
    }
  catch(const std::runtime_error &e)
    {
      std::cerr << "CFL Condition Error: " << e.what() << std::endl;
      throw; // Re-throw the exception to terminate the simulation
    }
}

void
MainSimulation::evaluateCarlemanNumber()
{
  try
    {
      double R = sim::utils::calculateCarlemanConvergenceNumber(
        initialConditions.getF0(), initialConditions.getF1(),
        initialConditions.getF2(), initialConditions.getU0s(),
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
