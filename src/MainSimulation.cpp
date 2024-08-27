#include "MainSimulation.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include <Eigen/Dense>
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
  checkStabilityConditions();

  initialConditions.computeInitialConditions();
  initialConditions.computeForcingBoundaryConditions(); // Compute the forcing
                                                        // boundary conditions
}

Eigen::MatrixXd
MainSimulation::prepareCarlemanMatrix()
{
  std::vector<int> dNs = matrix::calculateBlockSizes(params.N_max, params.nx);
  Eigen::MatrixXd  carlemanMatrix = matrix::assembleCarlemanMatrix(
     dNs, params.N_max, params.nx, params.ode_deg, F0, F1, F2);
  return carlemanMatrix;
}

void
MainSimulation::run()
{
  std::cout << "Running simulation..." << std::endl;
  std::cout << params << std::endl;
  std::cout << discretization << std::endl;
  // std::cout << initialConditions << std::endl;

  F0 = convertToDenseEigen(initialConditions.getF0());
  F1 = convertToDenseEigen(initialConditions.getF1());
  F2 = convertToDenseEigen(initialConditions.getF2());

  /*std::cout << "F0: " << F0 << std::endl;
  std::cout << "F1: " << F1 << std::endl;
  std::cout << "F2: " << F2 << std::endl;*/

  Eigen::EigenSolver<Eigen::MatrixXd> solver(F1);
  Eigen::VectorXcd                    eigenvalues = solver.eigenvalues();

  std::cout << "The eigenvalues of F1 are:\n" << eigenvalues << std::endl;
  evaluateCarlemanNumber();

  // Prepare the Carleman matrix
  Eigen::MatrixXd carlemanMatrix = prepareCarlemanMatrix();
  // print the size
  std::cout << "Carleman matrix size is: " << carlemanMatrix.rows() << " x "
            << carlemanMatrix.cols() << std::endl;
  // print the sum, the mean, and the max of the matrix
  std::cout << "Sum of the Carleman matrix: " << carlemanMatrix.sum()
            << std::endl;
  std::cout << "Mean of the Carleman matrix: " << carlemanMatrix.mean()
            << std::endl;
  std::cout << "Max of the Carleman matrix: " << carlemanMatrix.maxCoeff()
            << std::endl;

  // print the first 10 rows and columns of the matrix
  std::cout << "First 10 rows and columns of the Carleman matrix:\n"
            << carlemanMatrix.block(0, 0, 10, 10) << std::endl;

  // Solve the Carleman system
  carlemanSolver.solveCarlemanSystem(carlemanMatrix);

  // Solve the other systems
  eulerSolver.solveEuler(F0, F1, F2);
  // ode45Solver.solve();
  // pdeSolver.solve();
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

Eigen::MatrixXd
MainSimulation::convertToDenseEigen(const std::vector<std::vector<double>> &vec)
{
  // Get the dimensions of the input vector
  int rows = vec.size();
  int cols = vec[0].size();

  // Create an Eigen matrix with the same dimensions
  Eigen::MatrixXd eigenMatrix(rows, cols);

  // Copy the data from the 2D vector to the Eigen matrix
  for(int i = 0; i < rows; ++i)
    {
      for(int j = 0; j < cols; ++j)
        {
          eigenMatrix(i, j) = vec[i][j];
        }
    }

  return eigenMatrix;
}

} // namespace sim
