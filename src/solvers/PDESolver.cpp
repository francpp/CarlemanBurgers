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

double
PDESolver::F0_fun(double t, double x) const
{
  double L0 = params.L0;
  double U0 = params.U0;
  return U0 * exp(-pow(x - L0 / 4.0, 2) / (2 * pow(L0 / 32.0, 2))) *
         cos(2 * M_PI * t);
}

double
PDESolver::pde(double x, double t, double u, double dudx) const
{
  double nu = params.nu;
  double beta = params.beta;
  return nu * dudx - pow(u, 2) / 2.0 + F0_fun(t, x);
}

double
PDESolver::U0_fun(double x) const
{
  double L0 = params.L0;
  double U0 = params.U0;
  return -U0 * sin(2 * M_PI * x / L0);
}

void
PDESolver::apply_boundary_conditions(Eigen::VectorXd &u, double t) const
{
  // Assuming Dirichlet boundary conditions
  u(0) = 0;            // Left boundary
  u(u.size() - 1) = 0; // Right boundary
}

void
PDESolver::solvePDE(Eigen::MatrixXd &us_pde_full)
{
  int                        nx_pde = us_pde_full.cols();
  int                        nt_pde = us_pde_full.rows();
  const std::vector<double> &ts_pde = discretization.getTsPde();
  const std::vector<double> &xs_pde = discretization.getXsPde();
  double                     dx = xs_pde[1] - xs_pde[0];
  double                     dt = ts_pde[1] - ts_pde[0];

  // Initial condition
  for(int i = 0; i < nx_pde; ++i)
    {
      us_pde_full(0, i) = U0_fun(xs_pde[i]);
    }

  // Time-stepping loop with refinement
  for(int n = 0; n < nt_pde - 1; ++n)
    {
      Eigen::VectorXd u = us_pde_full.row(n);
      Eigen::VectorXd dudx(nx_pde);

      // Compute dudx using central difference with higher order accuracy
      for(int i = 1; i < nx_pde - 1; ++i)
        {
          dudx(i) = (u(i + 1) - u(i - 1)) / (2 * dx);
        }

      // Apply boundary conditions before updating
      apply_boundary_conditions(u, ts_pde[n]);

      // Update solution using explicit finite difference method with refinement
      for(int i = 1; i < nx_pde - 1; ++i)
        {
          double convection = -u(i) * (u(i + 1) - u(i - 1)) / (2 * dx);
          double diffusion =
            params.nu * (u(i + 1) - 2 * u(i) + u(i - 1)) / (dx * dx);
          double source = F0_fun(ts_pde[n], xs_pde[i]);

          us_pde_full(n + 1, i) = u(i) + dt * (convection + diffusion + source);
        }

      // Apply boundary conditions after updating
      Eigen::VectorXd next_u = us_pde_full.row(n + 1);
      apply_boundary_conditions(next_u, ts_pde[n]);

      // Update the next row in the matrix
      us_pde_full.row(n + 1) = next_u;
    }
}

void
PDESolver::solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1, Eigen::MatrixXd &F2)
{
  std::cout << "Solving with PDE solver" << std::endl;

  int nx_pde = params.nx_pde;
  int nt_pde = params.nt_pde;

  // Step 1: Solve the PDE on a coarse grid
  Eigen::MatrixXd us_pde_full(nt_pde, nx_pde);
  solvePDE(us_pde_full);

  // Step 2: Interpolate in time to match the finer time grid
  Eigen::MatrixXd us_pde_interp_temp(params.nt, nx_pde);
  for(int i = 0; i < params.nt; ++i)
    {
      us_pde_interp_temp.row(i) = Solver::interp1(
        discretization.getTsPde(), us_pde_full, discretization.getTs()[i]);
    }

  // Step 3: Interpolate in space to match the finer space grid
  Eigen::MatrixXd us_pde_fine(params.nt, params.nx);
  for(int k = 0; k < params.nx; ++k)
    {
      us_pde_fine.col(k) = Solver::interp1(discretization.getXsPde(),
                                           us_pde_interp_temp.transpose(),
                                           discretization.getXs()[k]);
    }

  // Store the final interpolated solution
  us_pde = us_pde_fine;
}

const Eigen::MatrixXd &
PDESolver::getUsPDE() const
{
  return us_pde;
}

} // namespace sim::solvers
