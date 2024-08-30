#include "boost/tuple/tuple.hpp"
#include "plots.hpp"

namespace sim
{
namespace plots
{

  Plotter::Plotter(
    const params::SimulationParameters          &params,
    const discretization::Discretization        &discretization,
    const initial_conditions::InitialConditions &initialConditions,
    const solvers::CarlemanSolver               &carlemanSolver,
    const solvers::EulerSolver                  &eulerSolver,
    const solvers::PDESolver                    &pdeSolver,
    const solvers::ODE45Solver                  &ode45Solver,
    const error_analysis::ErrorAnalysis         &errorAnalysis)
    : params(params), discretization(discretization),
      initialConditions(initialConditions), carlemanSolver(carlemanSolver),
      eulerSolver(eulerSolver), pdeSolver(pdeSolver), ode45Solver(ode45Solver),
      errorAnalysis(errorAnalysis)
  {}

  // Initialize plot settings
  void
  Plotter::initialize()
  {
    gp << "set terminal pngcairo enhanced size 1200,800\n";
  }

  void
  Plotter::plotCarlemanSolution()
  {
    // Set output file and initialize the plot title and axis labels
    gp << "set output 'solution_plot.png'\n";
    gp << "set title 'Final Carleman and PDE Solutions'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'u'\n";

    // Plot the Carleman solutions for N=1, N=2, N=3, N=4 and the PDE solution
    gp << "plot '-' with lines title 'Carleman Solution N=1', "
          "'-' with lines title 'Carleman Solution N=2', "
          "'-' with lines title 'Carleman Solution N=3', "
          "'-' with lines title 'Carleman Solution N=4', "
          "'-' with lines title 'PDE Solution'\n";

    // Extract the x-values (spatial grid points)
    const std::vector<double> &xs = discretization.getXs();
    size_t                     nx = xs.size(); // Number of spatial points

    // Get the Carleman solutions for the final time step
    const auto &us_c_N = carlemanSolver.getUsCN(); // Carleman solutions vector
    size_t      nt = us_c_N[0].rows();
    size_t      i_plot = nt - 1; // Index for the last timestamp

    // Send the data for each Carleman solution (N=1,2,3,4) at the final time
    // step
    for(int n = 0; n < 4; ++n)
      {
        std::vector<std::pair<double, double>> data_points(nx);
        for(size_t i = 0; i < nx; ++i)
          {
            data_points[i] = std::make_pair(xs[i], us_c_N[n](i_plot, i));
          }
        gp.send1d(data_points); // Carleman solution for N=n+1 at t=last
      }

    // Get the PDE solution for the final time step
    const Eigen::MatrixXd                 &us_pde = pdeSolver.getUsPDE();
    std::vector<std::pair<double, double>> pde_data_points(nx);
    for(size_t i = 0; i < nx; ++i)
      {
        pde_data_points[i] = std::make_pair(xs[i], us_pde(i_plot, i));
      }
    gp.send1d(pde_data_points); // PDE solution at t=last

    gp << "unset output\n";
  }

} // namespace plots
} // namespace sim
