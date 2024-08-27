#include "MainSimulation.hpp"
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
    errorAnalysis(eulerSolver, ode45Solver, pdeSolver, carlemanSolver, params)
{}

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
  F0 = matrixUtils::convertToDenseEigen(initialConditions.getF0());
  F1 = matrixUtils::convertToDenseEigen(initialConditions.getF1());
  F2 = matrixUtils::convertToDenseEigen(initialConditions.getF2());
}

void
MainSimulation::run()
{
  std::cout << "Running simulation..." << std::endl;
  std::cout << params << std::endl;
  std::cout << discretization << std::endl;
  // std::cout << initialConditions << std::endl;

  evaluateCarlemanNumber();

  // Solve the Carleman system
  carlemanSolver.solve(F0, F1, F2);

  // Solve the other systems
  eulerSolver.solve(F0, F1, F2);
  ode45Solver.solve(F0, F1, F2);
  pdeSolver.solve(F0, F1, F2);

  std::vector<Eigen::MatrixXd> us_c_N = carlemanSolver.getUsCN();
  Eigen::MatrixXd us_e = eulerSolver.getUsE();
  Eigen::MatrixXd us_d = ode45Solver.getUsD();
  Eigen::MatrixXd              us_pde = pdeSolver.getUsPDE();

  errorAnalysis.computeErrors();
  std::cout << "\nCD error: \n" << errorAnalysis.getEpsCDError() << std::endl;
  std::cout << "\nRelative CD error: \n"
            << errorAnalysis.getEpsRelCDError() << std::endl;
  std::cout << "\nCPDE error: \n"
            << errorAnalysis.getEpsCPDEError() << std::endl;
  std::cout << "\nRelative CPDE error: \n"
            << errorAnalysis.getEpsRelCPDEError() << std::endl;
  std::cout << "\nDPDE error: \n"
            << errorAnalysis.getEpsDPDEError() << std::endl;
  std::cout << "\nRelative DPDE error: \n"
            << errorAnalysis.getEpsRelDPDEError() << std::endl;
  std::cout << "\nDE error: \n" << errorAnalysis.getEpsDEError() << std::endl;
  std::cout << "\nDE error: \n"
            << errorAnalysis.getEpsRelDEError() << std::endl;
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

} // namespace sim
