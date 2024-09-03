#ifndef PDE_SOLVER_HPP
#define PDE_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include "utils/muparser_fun.hpp"
#include <Eigen/Dense>
#include <vector>

namespace sim::solvers
{
/**
 * @class PDESolver
 * @brief Implements a solver for partial differential equations (PDEs).
 *
 * The PDESolver class provides methods to solve PDEs using finite difference
 * methods with refinement, specifically designed for the simulation parameters
 * provided.
 */
class PDESolver : public Solver
{
public:
  /**
   * @brief Constructs a PDESolver with the given parameters, discretization,
   * and initial conditions.
   * @param params Reference to the simulation parameters.
   * @param discretization Reference to the discretization scheme.
   * @param initialConditions Reference to the initial conditions.
   */
  PDESolver(const params::SimulationParameters          &params,
            const discretization::Discretization        &discretization,
            const initial_conditions::InitialConditions &initialConditions);

  /**
   * @brief Solves the PDE system.
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
  const Eigen::MatrixXd &getUsPDE() const;

private:
  const params::SimulationParameters          &params;
  const discretization::Discretization        &discretization;
  const initial_conditions::InitialConditions &initialConditions;

  Eigen::MatrixXd us_pde; ///< Matrix to hold the solution for each time step.

  /**
   * @brief Forcing function \( F_0(t, x) \).
   * @param t Time variable.
   * @param x Spatial variable.
   * @return The value of the forcing function at time t and position x.
   */
  double F0_fun(double t, double x) const;

  /**
   * @brief The PDE function describing the system.
   * @param x Spatial variable.
   * @param t Time variable.
   * @param u Solution variable at position x and time t.
   * @param dudx Derivative of the solution with respect to x.
   * @return The value of the PDE at position x and time t.
   */
  double pde(double x, double t, double u, double dudx) const;

  /**
   * @brief Initial condition function \( U_0(x) \).
   * @param x Spatial variable.
   * @return The initial condition at position x.
   */
  double U0_fun(double x) const;

  /**
   * @brief Applies boundary conditions to the solution vector.
   * @param u Solution vector to which the boundary conditions are applied.
   * @param t Current time step.
   */
  void apply_boundary_conditions(Eigen::VectorXd &u, double t) const;

  /**
   * @brief Solves the PDE using finite difference methods.
   * @param us_pde_full Matrix to hold the full solution of the PDE.
   */
  void solvePDE(Eigen::MatrixXd &us_pde_full);
};

} // namespace sim::solvers

#endif // PDE_SOLVER_HPP
