#include "EulerSolver.hpp"
#include "matrix/kronp.hpp"
#include <iostream>

namespace sim::solvers
{

EulerSolver::EulerSolver(
  const params::SimulationParameters               &params,
  const sim::discretization::Discretization        &discretization,
  const sim::initial_conditions::InitialConditions &initialConditions)
  : params(params), discretization(discretization),
    initialConditions(initialConditions),
    us_e(Eigen::MatrixXd::Zero(params.nt, params.nx))
{}

void
EulerSolver::solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                   Eigen::MatrixXd &F2)
{
  double dt = discretization.getTs()[1] - discretization.getTs()[0];

  auto F0_interp = [&](double t) {
    Eigen::VectorXd interp = interp1(discretization.getTs(), F0, t);
    return interp;
  };

  auto burgers_odefun = [&](double t, const Eigen::MatrixXd &u) {
    Eigen::VectorXd burger_fun =
      F0_interp(t) + F1 * u + F2 * matrix::kron(u.sparseView(), u.sparseView());
    return burger_fun;
  };

  std::cout << "Solving direct Euler" << std::endl;

  Eigen::MatrixXd u0s = Eigen::Map<const Eigen::MatrixXd>(
    initialConditions.getU0s().data(), initialConditions.getU0s().size(), 1);

  us_e.row(0) = u0s.transpose();
  std::vector<double> ts = discretization.getTs();

  for(int k = 0; k < params.nt - 1; ++k)
    {
      us_e.row(k + 1) =
        us_e.row(k) +
        dt * burgers_odefun(ts[k], us_e.block(k, 0, 1, params.nx).transpose())
               .transpose()
               .row(0);
    }
}

const Eigen::MatrixXd &
EulerSolver::getUsE() const
{
  return us_e;
}

} // namespace sim::solvers
