#ifndef PLOTS_HPP
#define PLOTS_HPP

#include "error_analysis/ErrorAnalysis.hpp"
#include "gnuplot-iostream.hpp"
#include "solvers/CarlemanSolver.hpp"
#include "solvers/EulerSolver.hpp"
#include "solvers/ODE45Solver.hpp"
#include "solvers/PDESolver.hpp"
#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace sim::plots
{
class Plotter
{
public:
  Plotter(
    const sim::discretization::Discretization &discretization,
    const sim::error_analysis::ErrorAnalysis  &errorAnalysis,
    const sim::solvers::CarlemanSolver        &carlemanSolver,
    const sim::solvers::EulerSolver           &eulerSolver,
    const sim::solvers::ODE45Solver           &odeSolver,
    const sim::solvers::PDESolver             &pdeSolver,
    std::function<Eigen::MatrixXd(double, const Eigen::VectorXd &)> F0_fun);

  void generatePlots();

private:
  void plotInitialCondition();
  void plotErrorVsTime();
  void plotErrorConvergence();
  void formatPlot(Gnuplot &gp);

  int findIndex(const std::vector<double> &array, double value);
  const sim::discretization::Discretization &discretization;
  const sim::error_analysis::ErrorAnalysis  &errorAnalysis;
  const sim::solvers::CarlemanSolver        &carlemanSolver;
  const sim::solvers::EulerSolver           &eulerSolver;
  const sim::solvers::ODE45Solver           &odeSolver;
  const sim::solvers::PDESolver             &pdeSolver;
  std::function<Eigen::MatrixXd(double, const Eigen::VectorXd &)> F0_fun;
};

} // namespace sim::plots

#endif // PLOTS_HPP
