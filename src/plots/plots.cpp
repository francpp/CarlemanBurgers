#include "plots.hpp"
#include <Eigen/Dense>

namespace sim::plots
{

Plotter::Plotter(
  const sim::discretization::Discretization &discretization,
  const sim::error_analysis::ErrorAnalysis  &errorAnalysis,
  const sim::solvers::CarlemanSolver        &carlemanSolver,
  const sim::solvers::EulerSolver           &eulerSolver,
  const sim::solvers::ODE45Solver           &odeSolver,
  const sim::solvers::PDESolver             &pdeSolver,
  std::function<Eigen::MatrixXd(double, const Eigen::VectorXd &)> F0_fun)
  : discretization(discretization), errorAnalysis(errorAnalysis),
    carlemanSolver(carlemanSolver), eulerSolver(eulerSolver),
    odeSolver(odeSolver), pdeSolver(pdeSolver), F0_fun(F0_fun)
{}

void
Plotter::generatePlots()
{
  plotInitialCondition();
  plotErrorVsTime();
  plotErrorConvergence();
}

void
Plotter::plotInitialCondition()
{
  Gnuplot gp;

  auto us_pde = pdeSolver.getUsPDE();
  auto us_d = eulerSolver.getUsE();
  auto ts = discretization.getTs();
  auto xs = discretization.getXs();

  double t_plot = 1.0; // Example value, adjust as needed
  int    i_plot = findIndex(ts, t_plot);
  int    nx = xs.size();

  gp << "set multiplot layout 2,2 title 'Initial Condition and Solutions'\n";

  // Plot initial condition and source shape
  gp << "set title 'Initial condition and source shape'\n";
  gp << "plot '-' with lines title 'Initial condition', "
        "'-' with lines title 'Source shape', "
        "'-' with linespoints title 'Direct Euler solution'\n";
  gp.send1d(boost::make_tuple(xs, std::vector<double>(us_pde.row(0))));
  gp.send1d(boost::make_tuple(
    xs, std::vector<double>(
          F0_fun(1.0, Eigen::Map<Eigen::VectorXd>(xs.data(), xs.size())))));
  gp.send1d(boost::make_tuple(xs, std::vector<double>(us_d.row(i_plot))));

  // Plot Carleman solutions
  auto us_c_N = carlemanSolver.getUsCN();
  for(int N = 0; N < us_c_N.size(); ++N)
    {
      gp << "plot '-' with linespoints title 'Carleman solution, N=" << N
         << "'\n";
      gp.send1d(boost::make_tuple(xs, us_c_N[N].row(i_plot)));
    }
  gp << "unset multiplot\n";
}

void
Plotter::plotErrorVsTime()
{
  Gnuplot gp;
  auto    eps_c_d_N = errorAnalysis.getEpsCDError();
  auto    ts = discretization.getTs();

  gp << "set title 'Absolute error vs time'\n";
  gp << "set logscale y\n";
  gp << "plot ";

  for(int N = 0; N < eps_c_d_N.rows(); ++N)
    {
      if(N > 0)
        gp << ", ";
      gp << "'-' with lines title 'Carleman, N=" << N << "'";
    }
  gp << "\n";

  for(int N = 0; N < eps_c_d_N.rows(); ++N)
    {
      gp.send1d(boost::make_tuple(ts, eps_c_d_N.row(N)));
    }
}

void
Plotter::plotErrorConvergence()
{
  Gnuplot         gp;
  auto            eps_c_d_N = errorAnalysis.getEpsCDError();
  Eigen::VectorXd max_eps_c_d_N = eps_c_d_N.rowwise().maxCoeff();
  int             N_max = max_eps_c_d_N.size();

  gp << "set title 'Error Convergence'\n";
  gp << "set logscale y\n";
  gp << "plot '-' with linespoints title 'Time-maximum error'\n";
  gp.send1d(boost::make_tuple(Eigen::VectorXd::LinSpaced(N_max, 1, N_max),
                              max_eps_c_d_N));
}

void
Plotter::formatPlot(Gnuplot &gp)
{
  gp << "set grid\n";
  gp << "set xlabel 'x'\n";
  gp << "set ylabel 'u'\n";
}

int
Plotter::findIndex(const std::vector<double> &array, double value)
{
  for(int i = 0; i < array.size(); ++i)
    {
      if(array[i] >= value)
        return i;
    }
  return -1; // If not found, return an invalid index
}

} // namespace sim::plots
