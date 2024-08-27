#include "EulerSolver.hpp"
#include "MainSimulation.hpp"
#include "matrix/kronp.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace sim::solvers
{

EulerSolver::EulerSolver(
  const params::SimulationParameters               &params,
  const sim::discretization::Discretization        &discretization,
  const sim::initial_conditions::InitialConditions &initialConditions)
  : params(params), discretization(discretization),
    initialConditions(initialConditions)
{}

void
EulerSolver::solveEuler(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                        Eigen::MatrixXd &F2)
{
  double dt = discretization.getTs()[1] - discretization.getTs()[0];
  // Interpolator for F0
  auto F0_interp = [&](double t) {
    Eigen::VectorXd F0_interpolated = interp1(discretization.getTs(), F0, t);
    return F0_interpolated;
  };

  // Define the burgers_odefun lambda
  auto burgers_odefun = [&](double t, const Eigen::MatrixXd &u) {
    Eigen::MatrixXd output = F0_interp(t) + F1 * u + F2 * matrix::kron(u, u);
    return output;
  };

  std::cout << "Solving direct Euler" << std::endl;

  // Initialize the solution matrix
  Eigen::MatrixXd us_e = Eigen::MatrixXd::Zero(params.nt, params.nx);
  Eigen::MatrixXd u0s = Eigen::Map<const Eigen::MatrixXd>(
    initialConditions.getU0s().data(), initialConditions.getU0s().size(), 1);
  // Set the initial condition
  us_e.row(0) = u0s.transpose();
  std::vector<double> ts = discretization.getTs();

  // Time-stepping loop
  for(int k = 0; k < params.nt - 1; ++k)
    {
      us_e.row(k + 1) =
        us_e.row(k) +
        dt * burgers_odefun(ts[k], us_e.block(k, 0, 1, params.nx).transpose())
               .transpose()
               .row(0);
    }
  std::cout << us_e << std::endl;
}

Eigen::VectorXd
EulerSolver::interp1(const std::vector<double> &ts, const Eigen::MatrixXd &F0,
                     double t) const
{
  int n = ts.size();

  // Check if t is outside the range of ts
  if(t < ts.front() || t > ts.back())
    {
      throw std::out_of_range("t is out of bounds");
    }

  // Find the interval [t_i, t_i+1] such that t_i <= t < t_i+1
  int i = 0;
  while(i < n - 1 && t >= ts[i + 1])
    {
      ++i;
    }

  // Handle edge cases where t is exactly at the boundary
  if(t == ts[i])
    {
      return F0.row(i);
    }
  else if(t == ts[i + 1])
    {
      return F0.row(i + 1);
    }

  // Linear interpolation
  double t1 = ts[i];
  double t2 = ts[i + 1];
  double ratio = (t - t1) / (t2 - t1);

  return (1 - ratio) * F0.row(i) + ratio * F0.row(i + 1);
}
} // namespace sim::solvers
