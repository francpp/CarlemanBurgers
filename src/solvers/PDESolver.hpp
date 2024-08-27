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
class PDESolver
{
public:
  PDESolver(const params::SimulationParameters          &params,
            const discretization::Discretization        &discretization,
            const initial_conditions::InitialConditions &initialConditions);

  void solvePDE(Eigen::MatrixXd &F0);

  const Eigen::MatrixXd &
  getUsPDE() const; // Method to access the solution matrix

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_pde; // Member to hold the solution matrix

  Eigen::VectorXd interp1(const std::vector<double> &ts,
                          const Eigen::MatrixXd &F0, double t) const;
  double          F0Function(double t, double x);
};
} // namespace sim::solvers

#endif // PDE_SOLVER_HPP
