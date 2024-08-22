#include "MainSimulation.hpp"

namespace sim
{

MainSimulation::MainSimulation(params::SimulationParameters &params)
  : params(params), // Initialize the member variable params with the reference
    discretization(params), initialConditions(params, discretization),
    matrixOperations(params, discretization),
    eulerSolver(params, discretization, initialConditions, matrixOperations),
    ode45Solver(params, discretization, initialConditions, matrixOperations),
    pdeSolver(params, discretization, initialConditions, matrixOperations),
    carlemanSolver(params, discretization, initialConditions, matrixOperations),
    errorAnalysis(eulerSolver, ode45Solver, pdeSolver, carlemanSolver)
{
  // Constructor implementation
}

void
MainSimulation::initialize()
{
  params
    .initialize(); // Initialize params (which will now affect all components)
}

void
MainSimulation::run()
{
  std::cout << "Running simulation..." << std::endl;
  std::cout << params << std::endl;
  discretization.createDiscretization();
  std::cout << discretization << std::endl;

  initialConditions.computeInitialConditions();

  checkStabilityConditions();
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

} // namespace sim
