#include "ErrorAnalysis.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

namespace sim::error_analysis
{

ErrorAnalysis::ErrorAnalysis(const sim::solvers::EulerSolver    &eulerSolver,
                             const sim::solvers::ODE45Solver    &odeSolver,
                             const sim::solvers::PDESolver      &pdeSolver,
                             const sim::solvers::CarlemanSolver &carlemanSolver,
                             const params::SimulationParameters &params)
  : eulerSolver(eulerSolver), odeSolver(odeSolver), pdeSolver(pdeSolver),
    carlemanSolver(carlemanSolver), params(params),
    eps_c_d_N(Eigen::MatrixXd::Zero(params.N_max, params.nt)),
    eps_rel_c_d_N(Eigen::MatrixXd::Zero(params.N_max, params.nt)),
    eps_c_pde_N(Eigen::MatrixXd::Zero(params.N_max, params.nt)),
    eps_rel_c_pde_N(Eigen::MatrixXd::Zero(params.N_max, params.nt)),
    eps_d_pde(Eigen::VectorXd::Zero(params.nt)),
    eps_rel_d_pde(Eigen::VectorXd::Zero(params.nt)),
    eps_d_e(Eigen::VectorXd::Zero(params.nt)),
    eps_rel_d_e(Eigen::VectorXd::Zero(params.nt))
{}

void
ErrorAnalysis::computeErrors()
{
  auto us_c_N = carlemanSolver.getUsCN();
  auto us_e = eulerSolver.getUsE();
  auto us_d = odeSolver.getUsD();
  auto us_pde = pdeSolver.getUsPDE();

  for(int N = 0; N < params.N_max; ++N)
    {
      // Absolute error between Carleman solution and ODE solution
      Eigen::MatrixXd dus_c_d_N = us_c_N[N] - us_d;

      // Relative error between Carleman solution and ODE solution
      Eigen::MatrixXd dus_rel_c_d_N = dus_c_d_N.array() / us_d.array();
      dus_rel_c_d_N = dus_rel_c_d_N.unaryExpr(
        [](double v) { return std::isfinite(v) ? v : 0; });

      // Absolute error between Carleman solution and PDE solution
      Eigen::MatrixXd dus_c_pde_N = us_c_N[N] - us_pde;

      // Relative error between Carleman solution and PDE solution
      Eigen::MatrixXd dus_rel_c_pde_N = dus_c_pde_N.array() / us_pde.array();
      dus_rel_c_pde_N = dus_rel_c_pde_N.unaryExpr(
        [](double v) { return std::isfinite(v) ? v : 0; });

      // Compute norms over space for each time step
      for(int k = 0; k < params.nt; ++k)
        {
          // L2 norm of the absolute error between Carleman and ODE solutions
          eps_c_d_N(N, k) = dus_c_d_N.row(k).norm();

          // L-inf norm of the relative error between Carleman and ODE solutions
          eps_rel_c_d_N(N, k) = dus_rel_c_d_N.row(k).lpNorm<Eigen::Infinity>();

          // L2 norm of the absolute error between Carleman and PDE solutions
          eps_c_pde_N(N, k) = dus_c_pde_N.row(k).norm();

          // L-inf norm of the relative error between Carleman and PDE solutions
          eps_rel_c_pde_N(N, k) =
            dus_rel_c_pde_N.row(k).lpNorm<Eigen::Infinity>();
        }
    }

  // Absolute error between ODE and PDE solutions
  Eigen::MatrixXd dus_d_pde = us_d - us_pde;

  // Relative error between ODE and PDE solutions
  Eigen::MatrixXd dus_rel_d_pde = dus_d_pde.array() / us_pde.array();
  dus_rel_d_pde =
    dus_rel_d_pde.unaryExpr([](double v) { return std::isfinite(v) ? v : 0; });

  // Absolute error between ODE and Euler solutions
  Eigen::MatrixXd dus_d_e = us_d - us_e;

  // Relative error between ODE and Euler solutions
  Eigen::MatrixXd dus_rel_d_e = dus_d_e.array() / us_e.array();
  dus_rel_d_e =
    dus_rel_d_e.unaryExpr([](double v) { return std::isfinite(v) ? v : 0; });

  // Compute norms over space for each time step
  for(int k = 0; k < params.nt; ++k)
    {
      // L2 norm of the absolute error between ODE and PDE solutions
      eps_d_pde(k) = dus_d_pde.row(k).norm();

      // L-inf norm of the relative error between ODE and PDE solutions
      eps_rel_d_pde(k) = dus_rel_d_pde.row(k).lpNorm<Eigen::Infinity>();

      // L2 norm of the absolute error between ODE and Euler solutions
      eps_d_e(k) = dus_d_e.row(k).norm();

      // L-inf norm of the relative error between ODE and Euler solutions
      eps_rel_d_e(k) = dus_rel_d_e.row(k).lpNorm<Eigen::Infinity>();
    }
}

// Getter methods to retrieve computed errors

const Eigen::MatrixXd &
ErrorAnalysis::getEpsCDError() const
{
  return eps_c_d_N;
}

const Eigen::MatrixXd &
ErrorAnalysis::getEpsRelCDError() const
{
  return eps_rel_c_d_N;
}

const Eigen::MatrixXd &
ErrorAnalysis::getEpsCPDEError() const
{
  return eps_c_pde_N;
}

const Eigen::MatrixXd &
ErrorAnalysis::getEpsRelCPDEError() const
{
  return eps_rel_c_pde_N;
}

const Eigen::VectorXd &
ErrorAnalysis::getEpsDPDEError() const
{
  return eps_d_pde;
}

const Eigen::VectorXd &
ErrorAnalysis::getEpsRelDPDEError() const
{
  return eps_rel_d_pde;
}

const Eigen::VectorXd &
ErrorAnalysis::getEpsDEError() const
{
  return eps_d_e;
}

const Eigen::VectorXd &
ErrorAnalysis::getEpsRelDEError() const
{
  return eps_rel_d_e;
}

} // namespace sim::error_analysis
