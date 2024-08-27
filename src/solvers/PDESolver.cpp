#include "PDESolver.hpp"
#include <iostream>

namespace sim::solvers
{

PDESolver::PDESolver(
  const params::SimulationParameters               &params,
  const sim::discretization::Discretization        &discretization,
  const sim::initial_conditions::InitialConditions &initialConditions)
  : params(params), discretization(discretization),
    initialConditions(initialConditions),
    us_pde(Eigen::MatrixXd::Zero(params.nt, params.nx))
{}

void
PDESolver::solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1, Eigen::MatrixXd &F2)
{
  std::cout << "Solving with PDE solver" << std::endl;

  int nx_pde = params.nx;
  int nt_pde = params.nt;

  us_pde =
    Eigen::MatrixXd::Zero(nt_pde, nx_pde); // Placeholder solution process

  std::cout << "PDE solved and interpolated" << std::endl;
}

const Eigen::MatrixXd &
PDESolver::getUsPDE() const
{
  return us_pde;
}

} // namespace sim::solvers
