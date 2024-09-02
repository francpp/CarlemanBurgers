#ifndef PLOTS_HPP
#define PLOTS_HPP

#include "discretization/Discretization.hpp"
#include "error_analysis/ErrorAnalysis.hpp"
#include "gnuplot-iostream.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include "solvers/CarlemanSolver.hpp"
#include "solvers/EulerSolver.hpp"
#include "solvers/ODE45Solver.hpp"
#include "solvers/PDESolver.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace sim
{
namespace plots
{

  class Plotter
  {
  public:
    // Constructor
    Plotter(const params::SimulationParameters          &params,
            const discretization::Discretization        &discretization,
            const initial_conditions::InitialConditions &initialConditions,
            const solvers::CarlemanSolver               &carlemanSolver,
            const solvers::EulerSolver                  &eulerSolver,
            const solvers::PDESolver                    &pdeSolver,
            const solvers::ODE45Solver                  &ode45Solver,
            const error_analysis::ErrorAnalysis         &errorAnalysis);

    // Method to initialize plot settings
    void initialize();

    // Method to plot solutions
    void plotSolution();

    // Method to plot the absolute error
    void plotErrors();

    // Method to plot the error convergence
    void plotErrorConvergence();

  private:
    // Member variables for plot settings and parameters
    Gnuplot gp;
    const params::SimulationParameters
      &params; // Reference to the simulation parameters
    const discretization::Discretization
      &discretization; // Reference to the discretization
    const initial_conditions::InitialConditions
      &initialConditions; // Reference to the initial conditions
    const solvers::CarlemanSolver
      &carlemanSolver;                       // Reference to the Carleman solver
    const solvers::EulerSolver &eulerSolver; // Reference to the Euler solver
    const solvers::PDESolver   &pdeSolver;   // Reference to the PDE solver
    const solvers::ODE45Solver &ode45Solver; // Reference to the ODE45 solver
    const error_analysis::ErrorAnalysis
      &errorAnalysis; // Reference to the error analysis
  };

} // namespace plots
} // namespace sim
#endif // PLOTS_HPP
