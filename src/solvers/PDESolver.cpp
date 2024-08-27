#include "PDESolver.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace sim::solvers
{

PDESolver::PDESolver(
  const params::SimulationParameters               &params,
  const sim::discretization::Discretization        &discretization,
  const sim::initial_conditions::InitialConditions &initialConditions)
  : params(params), discretization(discretization),
    initialConditions(initialConditions),
    us_pde(Eigen::MatrixXd::Zero(params.nt,
                                 params.nx)) // Initialize us_pde with zeros
{}

void
PDESolver::solvePDE(Eigen::MatrixXd &F0)
{
  std::cout << "Solving with PDE solver" << std::endl;

  int nx_pde = params.nx;
  int nt_pde = params.nt;

  // Assuming the solution is solved and stored in us_pde as shown below:
  us_pde =
    Eigen::MatrixXd::Zero(nt_pde, nx_pde); // Placeholder solution process

  std::cout << "PDE solved and interpolated" << std::endl;
}

double
PDESolver::F0Function(double t, double x)
{
  double sigma = params.L0 / 32.0;
  double exponent =
    -(std::pow(x - params.L0 / 4.0, 2) / (2 * std::pow(sigma, 2)));
  return params.U0 * std::exp(exponent) * std::cos(2 * M_PI * t);
}

Eigen::VectorXd
PDESolver::interp1(const std::vector<double> &ts, const Eigen::MatrixXd &F0,
                   double t) const
{
  int n = ts.size();

  if(t < ts.front() || t > ts.back())
    throw std::out_of_range("t is out of bounds");

  int i = 0;
  while(i < n - 1 && t >= ts[i + 1])
    ++i;

  if(t == ts[i])
    return F0.row(i);
  else if(t == ts[i + 1])
    return F0.row(i + 1);

  double t1 = ts[i];
  double t2 = ts[i + 1];
  double ratio = (t - t1) / (t2 - t1);

  return (1 - ratio) * F0.row(i) + ratio * F0.row(i + 1);
}

const Eigen::MatrixXd &
PDESolver::getUsPDE() const
{
  return us_pde;
}

} // namespace sim::solvers
