#include "boost/tuple/tuple.hpp"
#include "plots.hpp"
#include <boost/filesystem.hpp>
#include <cmath>

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

  void
  Plotter::initialize()
  {
    // Create the output directory if it doesn't exist
    boost::filesystem::path output_dir("output");
    if(!boost::filesystem::exists(output_dir))
      {
        boost::filesystem::create_directory(output_dir);
      }

    gp << "set terminal pngcairo enhanced size 1200,800\n";
  }

  void
  Plotter::plotSolution()
  {
    // Set output file and initialize the plot title and axis labels
    gp << "set output 'output/solution_plot.png'\n";
    gp << "set title 'Solution at T_{nl}/3 and Initial Condition'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'u'\n";
    gp << "set key left top\n"; // Move legend to top left for better visibility

    // Plot the Carleman solutions, initial condition, and PDE solution
    gp << "plot '-' with lines title 'Initial Condition', "
          "'-' with lines title 'Source Shape', "
          "'-' with linespoints title 'Direct Euler at T_{nl}/3', "
          "'-' with lines title 'Carleman N=1 at T_{nl}/3', "
          "'-' with lines title 'Carleman N=N_{max} at T_{nl}/3'\n";

    const std::vector<double> &xs = discretization.getXs();
    size_t                     nx = xs.size();

    // Plot the initial condition
    std::vector<std::pair<double, double>> init_condition(nx);
    for(size_t i = 0; i < nx; ++i)
      {
        init_condition[i] = std::make_pair(xs[i], 0);
      }
    gp.send1d(init_condition);

    // Plot the source shape
    std::vector<std::pair<double, double>> source_shape(nx);
    for(size_t i = 0; i < nx; ++i)
      {
        source_shape[i] = std::make_pair(xs[i], 0);
      }
    gp.send1d(source_shape);

    // Plot the Direct Euler solution at T_{nl}/3
    double              T_nl = params.L0 / params.U0;
    double              t_plot = T_nl / 3;
    std::vector<double> ts = discretization.getTs();
    int                 i_plot = -1;

    // Find i_plot
    for(size_t i = 0; i < ts.size(); ++i)
      {
        if(ts[i] >= t_plot)
          {
            i_plot = i + 1;
            break;
          }
      }

    std::vector<std::pair<double, double>> euler_solution(nx);
    const Eigen::MatrixXd &us_d = ode45Solver.getUsD(); // ODE45 solution
    for(size_t i = 0; i < nx; ++i)
      {
        euler_solution[i] = std::make_pair(xs[i], us_d(i_plot, i));
      }
    gp.send1d(euler_solution);

    // Plot Carleman solutions for N=1 and N=N_max at T_{nl}/3
    const auto &us_c_N = carlemanSolver.getUsCN();
    for(int n : {0, params.N_max - 1})
      {
        std::vector<std::pair<double, double>> carleman_solution(nx);
        for(size_t i = 0; i < nx; ++i)
          {
            carleman_solution[i] = std::make_pair(xs[i], us_c_N[n](i_plot, i));
          }
        gp.send1d(carleman_solution);
      }

    gp << "unset output\n";
  }

  void
  Plotter::plotErrors()
  {
    gp << "set output 'output/error_plot.png'\n";
    gp << "set title 'Absolute L_2 Error between Carleman and ODE45 "
          "Solutions'\n";
    gp << "set xlabel 't'\n";
    gp << "set ylabel '||\\varepsilon_{\\mathrm{abs}}||_2'\n";
    gp << "set logscale y\n"; // Use logarithmic scale on y-axis
    gp << "plot ";

    // Define the time vector
    const std::vector<double> &ts = discretization.getTs();
    size_t                     nt = ts.size();
    double                     T_nl = params.L0 / params.U0;
    double                     t_plot = T_nl / 3;
    int                        i_plot = -1;

    // Find i_plot
    for(size_t i = 0; i < ts.size(); ++i)
      {
        if(ts[i] >= t_plot)
          {
            i_plot = i;
            break;
          }
      }

    // Prepare Carleman errors for different N
    const auto &eps_c_d_N = errorAnalysis.getEpsCDError();
    for(size_t N = 0; N < eps_c_d_N.rows(); ++N)
      {
        if(N > 0)
          gp << ", "; // Separate multiple plots
        gp << "'-' with lines title 'Carleman N=" << (N + 1) << "'";
      }
    gp << ", '-' with lines lt -1 title 'T_{nl}/3'\n";

    // Send the data for the errors
    for(size_t i = 0; i < eps_c_d_N.rows(); ++i)
      {
        std::vector<std::pair<double, double>> error_data(nt);
        for(size_t k = 0; k < nt; ++k)
          {
            error_data[k] = std::make_pair(ts[k], eps_c_d_N(i, k));
          }
        gp.send1d(error_data);
      }

    // Add the vertical line at T_{nl}/3
    double min_val =
      eps_c_d_N.block(0, i_plot, eps_c_d_N.rows(), eps_c_d_N.cols() - i_plot)
        .minCoeff();
    double max_val = eps_c_d_N.row(0).maxCoeff();
    std::vector<std::pair<double, double>> t_plot_line = {
      {t_plot, 0.99 * min_val}, {t_plot, 1.1 * max_val}};
    gp.send1d(t_plot_line);

    gp << "unset output\n";
  }

  void
  Plotter::plotErrorConvergence()
  {
    gp << "set output 'output/error_convergence_plot.png'\n";
    gp << "set title 'Error Convergence of Carleman Solutions'\n";
    gp << "set xlabel 'N'\n";
    gp << "set ylabel 'max_t ||\\varepsilon_{\\mathrm{abs}}||_2'\n";
    gp << "set logscale y\n";

    const auto     &eps_c_d_N = errorAnalysis.getEpsCDError();
    size_t          N_max = eps_c_d_N.rows();
    Eigen::VectorXd max_values;
    max_values = eps_c_d_N.rowwise().maxCoeff();

    // Prepare data for max error over time for each N
    std::vector<std::pair<int, double>> max_error_data(N_max);
    for(size_t N = 0; N < N_max; ++N)
      {
        max_error_data[N] = std::make_pair(N + 1, max_values[N]);
      }

    // Plot the maximum error vs. N
    gp << "plot '-' with linespoints title 'Error Convergence'\n";
    gp.send1d(max_error_data);

    gp << "unset output\n";
  }

} // namespace plots
} // namespace sim
