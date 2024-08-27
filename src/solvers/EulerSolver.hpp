#ifndef EULER_SOLVER_HPP
#define EULER_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"

namespace sim::solvers
{
class EulerSolver : public Solver
{
public:
  EulerSolver(const params::SimulationParameters          &params,
              const discretization::Discretization        &discretization,
              const initial_conditions::InitialConditions &initialConditions);

  void solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
             Eigen::MatrixXd &F2) override;

  const Eigen::MatrixXd &getUsE() const;

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_e;
};
} // namespace sim::solvers

#endif // EULER_SOLVER_HPP
