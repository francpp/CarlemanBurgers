#ifndef EULER_SOLVER_HPP
#define EULER_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include <Eigen/Dense>
#include <vector>

namespace sim::solvers
{
class EulerSolver
{
public:
  EulerSolver(const params::SimulationParameters          &params,
              const discretization::Discretization        &discretization,
              const initial_conditions::InitialConditions &initiaConditions);

  void solveEuler(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                  Eigen::MatrixXd &F2);

  const Eigen::MatrixXd &getUsE() const; // Method to access the solution matrix

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_e; // Member to hold the solution matrix

  Eigen::VectorXd interp1(const std::vector<double> &ts,
                          const Eigen::MatrixXd &F0, double t) const;
};
} // namespace sim::solvers

#endif // EULER_SOLVER_HPP
