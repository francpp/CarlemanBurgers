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

  // initialConditions.computeInitialConditions();
  // std::cout << "Discretization and initial conditions set up." << std::endl;
}

} // namespace sim
