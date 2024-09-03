#ifndef ODE45_SOLVER_HPP
#define ODE45_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"

namespace sim::solvers
{
/**
 * @class ODE45Solver
 * @brief Implements the Runge-Kutta method (RK45) for solving differential
 * equations.
 *
 * The ODE45Solver class provides an implementation of the Runge-Kutta method of
 * order 4(5) (ODE45) to solve systems of differential equations, specifically
 * designed for the simulation parameters provided.
 */
class ODE45Solver : public Solver
{
public:
  /**
   * @brief Constructs an ODE45Solver with the given parameters, discretization,
   * and initial conditions.
   * @param params Reference to the simulation parameters.
   * @param discretization Reference to the discretization scheme.
   * @param initialConditions Reference to the initial conditions.
   */
  ODE45Solver(const params::SimulationParameters          &params,
              const discretization::Discretization        &discretization,
              const initial_conditions::InitialConditions &initialConditions);

  /**
   * @brief Solves the system using the Runge-Kutta method (ODE45).
   *
   * This function overrides the pure virtual solve function from the Solver
   * interface.
   *
   * @param F0 Forcing matrix F0.
   * @param F1 Linear coefficient matrix F1.
   * @param F2 Quadratic interaction matrix F2.
   */
  void solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
             Eigen::MatrixXd &F2) override;

  /**
   * @brief Returns the computed solution matrix.
   * @return A reference to the matrix containing the solution.
   */
  const Eigen::MatrixXd &getUsD() const;

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_d; ///< Matrix to hold the solution for each time step.
};

} // namespace sim::solvers

#endif // ODE45_SOLVER_HPP
