#include "ODE45Solver.hpp"
#include "matrix/kronp.hpp"
#include <iostream>

namespace sim::solvers
{

ODE45Solver::ODE45Solver(
  const params::SimulationParameters               &params,
  const sim::discretization::Discretization        &discretization,
  const sim::initial_conditions::InitialConditions &initialConditions)
  : params(params), discretization(discretization),
    initialConditions(initialConditions),
    us_d(Eigen::MatrixXd::Zero(params.nt, params.nx))
{}

void
ODE45Solver::solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                   Eigen::MatrixXd &F2)
{
  double dt_ode = discretization.getTsOde()[1] - discretization.getTsOde()[0];

  auto F0_interp = [&](double t) {
    Eigen::VectorXd interp = interp1(discretization.getTs(), F0, t);
    return interp;
  };

  auto burgers_odefun = [&](double t, const Eigen::MatrixXd &u) {
    Eigen::VectorXd burger_fun =
      F0_interp(t) + F1 * u + F2 * matrix::kron(u, u);
    return burger_fun;
  };

  std::cout << "Solving with Runge-Kutta method (RK45)" << std::endl;

  Eigen::MatrixXd us_ode = Eigen::MatrixXd::Zero(10 * params.nt, params.nx);
  Eigen::MatrixXd u0s = Eigen::Map<const Eigen::MatrixXd>(
    initialConditions.getU0s().data(), initialConditions.getU0s().size(), 1);

  us_ode.row(0) = u0s.transpose();
  std::vector<double> ts_ode = discretization.getTsOde();
  double              t = ts_ode[0];

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

  std::vector<double> ts = discretization.getTs();

  for(int i = 0; i < params.nt; ++i)
    {
      Eigen::VectorXd u_interpolated = interp1(ts_ode, us_ode, ts[i]);
      us_d.row(i) = u_interpolated;
    }
}

const Eigen::MatrixXd &
ODE45Solver::getUsD() const
{
  return us_d;
}

} // namespace sim::solvers
