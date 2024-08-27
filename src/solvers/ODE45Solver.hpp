#ifndef ODE45_SOLVER_HPP
#define ODE45_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"

namespace sim::solvers
{
class ODE45Solver : public Solver
{
public:
  ODE45Solver(const params::SimulationParameters          &params,
              const discretization::Discretization        &discretization,
              const initial_conditions::InitialConditions &initialConditions);

  void solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
             Eigen::MatrixXd &F2) override;

  const Eigen::MatrixXd &getUsD() const;

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_d;
};
} // namespace sim::solvers

#endif // ODE45_SOLVER_HPP
