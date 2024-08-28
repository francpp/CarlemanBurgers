#ifndef PDE_SOLVER_HPP
#define PDE_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include <Eigen/Dense>
#include <vector>

namespace sim::solvers
{
class PDESolver : public Solver
{
public:
  PDESolver(const params::SimulationParameters          &params,
            const discretization::Discretization        &discretization,
            const initial_conditions::InitialConditions &initialConditions);

  void solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
             Eigen::MatrixXd &F2) override;

  const Eigen::MatrixXd &getUsPDE() const;

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_pde;

  double F0_fun(double t, double x) const;
  double pde(double x, double t, double u, double dudx) const;
  double initial_condition(double x) const;
  void   apply_boundary_conditions(Eigen::VectorXd &u, double t) const;

  void solvePDE(Eigen::MatrixXd &us_pde_full);
};
} // namespace sim::solvers

#endif // PDE_SOLVER_HPP
