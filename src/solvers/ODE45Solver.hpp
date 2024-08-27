#ifndef ODE45_SOLVER_HPP
#define ODE45_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include <Eigen/Dense>
#include <vector>

namespace sim::solvers
{
class ODE45Solver
{
public:
  ODE45Solver(const params::SimulationParameters          &params,
              const discretization::Discretization        &discretization,
              const initial_conditions::InitialConditions &initiaConditions);

  void solveODE45(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                  Eigen::MatrixXd &F2);

  const Eigen::MatrixXd &getUsD() const; // Method to access the solution matrix

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_d; // Member to hold the solution matrix

  Eigen::VectorXd interp1(const std::vector<double> &ts,
                          const Eigen::MatrixXd &F0, double t) const;
};
} // namespace sim::solvers

#endif // ODE45_SOLVER_HPP
