#include "MainSimulation.hpp"
#include "ODE45Solver.hpp"
#include "matrix/kronp.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace sim::solvers
{

ODE45Solver::ODE45Solver(
  const params::SimulationParameters               &params,
  const sim::discretization::Discretization        &discretization,
  const sim::initial_conditions::InitialConditions &initialConditions)
  : params(params), discretization(discretization),
    initialConditions(initialConditions),
    us_d(
      Eigen::MatrixXd::Zero(params.nt, params.nx)) // Initialize us_d with zeros
{}

void
ODE45Solver::solveODE45(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                        Eigen::MatrixXd &F2)
{
  double dt_ode = discretization.getTsOde()[1] - discretization.getTsOde()[0];

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

  std::cout << "Solving with Runge-Kutta method (RK45)" << std::endl;

  // Initialize the solution matrix for the ODE solver
  Eigen::MatrixXd us_ode = Eigen::MatrixXd::Zero(10 * params.nt, params.nx);
  Eigen::MatrixXd u0s = Eigen::Map<const Eigen::MatrixXd>(
    initialConditions.getU0s().data(), initialConditions.getU0s().size(), 1);

  // Set the initial condition
  us_ode.row(0) = u0s.transpose();

  std::vector<double> ts_ode = discretization.getTsOde();
  double              t = ts_ode[0];

  // Time stepping loop using a simple Runge-Kutta 4th order (RK4)
  for(int i = 1; i < ts_ode.size(); ++i)
    {
      double h = ts_ode[i] - ts_ode[i - 1];

      Eigen::MatrixXd u = us_ode.row(i - 1).transpose();

      Eigen::MatrixXd k1 = burgers_odefun(t, u);
      Eigen::MatrixXd k2 = burgers_odefun(t + 0.5 * h, u + 0.5 * h * k1);
      Eigen::MatrixXd k3 = burgers_odefun(t + 0.5 * h, u + 0.5 * h * k2);
      Eigen::MatrixXd k4 = burgers_odefun(t + h, u + h * k3);

      Eigen::MatrixXd u_next = u + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
      us_ode.row(i) = u_next.transpose();
      t = ts_ode[i];
    }

  // At this point, us_ode contains the solution for each time step

  // Interpolating the ODE solution to the original time grid
  std::vector<double> ts = discretization.getTs();

  for(int i = 0; i < params.nt; ++i)
    {
      Eigen::VectorXd u_interpolated = interp1(ts_ode, us_ode, ts[i]);
      us_d.row(i) = u_interpolated;
    }
  std::cout << "ODE solved" << std::endl;
}

Eigen::VectorXd
ODE45Solver::interp1(const std::vector<double> &ts, const Eigen::MatrixXd &F0,
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

const Eigen::MatrixXd &
ODE45Solver::getUsD() const
{
  return us_d;
}

} // namespace sim::solvers
