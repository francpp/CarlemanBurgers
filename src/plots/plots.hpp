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
  /**
   * @class Plotter
   * @brief Provides functionality for plotting solutions, errors, and error
   * convergence.
   *
   * The Plotter class utilizes gnuplot to generate plots for solutions, errors,
   * and the convergence of errors in the context of the simulation.
   */
  class Plotter
  {
  public:
    /**
     * @brief Constructs a Plotter with the given simulation parameters and
     * solvers.
     * @param params Reference to the simulation parameters.
     * @param discretization Reference to the discretization scheme.
     * @param initialConditions Reference to the initial conditions.
     * @param carlemanSolver Reference to the Carleman solver.
     * @param eulerSolver Reference to the Euler solver.
     * @param pdeSolver Reference to the PDE solver.
     * @param ode45Solver Reference to the ODE45 solver.
     * @param errorAnalysis Reference to the error analysis.
     */
    Plotter(const params::SimulationParameters          &params,
            const discretization::Discretization        &discretization,
            const initial_conditions::InitialConditions &initialConditions,
            const solvers::CarlemanSolver               &carlemanSolver,
            const solvers::EulerSolver                  &eulerSolver,
            const solvers::PDESolver                    &pdeSolver,
            const solvers::ODE45Solver                  &ode45Solver,
            const error_analysis::ErrorAnalysis         &errorAnalysis);

    /**
     * @brief Initializes the plot settings and creates the output directory.
     */
    void initialize();

    /**
     * @brief Plots the solution at a specific time point.
     */
    void plotSolution();

    /**
     * @brief Plots the absolute error between the Carleman and ODE45 solutions.
     */
    void plotErrors();

    /**
     * @brief Plots the convergence of the error as the Carleman order
     * increases.
     */
    void plotErrorConvergence();

  private:
    Gnuplot gp; ///< Gnuplot object for plotting.

    const params::SimulationParameters
      &params; ///< Reference to the simulation parameters.
    const discretization::Discretization
      &discretization; ///< Reference to the discretization.
    const initial_conditions::InitialConditions
      &initialConditions; ///< Reference to the initial conditions.
    const solvers::CarlemanSolver
      &carlemanSolver; ///< Reference to the Carleman solver.
    const solvers::EulerSolver &eulerSolver; ///< Reference to the Euler solver.
    const solvers::PDESolver   &pdeSolver;   ///< Reference to the PDE solver.
    const solvers::ODE45Solver &ode45Solver; ///< Reference to the ODE45 solver.
    const error_analysis::ErrorAnalysis
      &errorAnalysis; ///< Reference to the error analysis.
  };

} // namespace plots
} // namespace sim

#endif // PLOTS_HPP
