#ifndef CARLEMAN_SOLVER_HPP
#define CARLEMAN_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace sim
{
namespace solvers
{
  /**
   * @class CarlemanSolver
   * @brief Solver for systems using Carleman linearization.
   *
   * The CarlemanSolver class implements the Solver interface and provides
   * methods to solve a system using Carleman linearization, which is used
   * to approximate solutions of nonlinear differential equations by
   * transforming them into a linear system of equations.
   */
  class CarlemanSolver : public Solver
  {
  public:
    /**
     * @brief Constructs a CarlemanSolver with the given parameters,
     * discretization, and initial conditions.
     * @param params Reference to the simulation parameters.
     * @param discretization Reference to the discretization scheme.
     * @param initialConditions Reference to the initial conditions.
     */
    CarlemanSolver(
      const params::SimulationParameters          &params,
      const discretization::Discretization        &discretization,
      const initial_conditions::InitialConditions &initialConditions);

    /**
     * @brief Solves the system using Carleman linearization.
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
     * @brief Returns the computed solution matrices.
     * @return A vector of matrices, each representing the solution for a
     * different Carleman order.
     */
    const std::vector<Eigen::MatrixXd> &getUsCN() const;

  private:
    /**
     * @brief Prepares the Carleman matrix for the solver.
     * @param F0 Forcing matrix F0.
     * @param F1 Linear coefficient matrix F1.
     * @param F2 Quadratic interaction matrix F2.
     * @return The assembled Carleman matrix as a sparse matrix.
     */
    Eigen::SparseMatrix<double> prepareCarlemanMatrix(Eigen::MatrixXd &F0,
                                                      Eigen::MatrixXd &F1,
                                                      Eigen::MatrixXd &F2);

    const params::SimulationParameters          &params;
    const discretization::Discretization        &discretization;
    const initial_conditions::InitialConditions &initialConditions;

    std::vector<Eigen::MatrixXd>
      us_c_N; ///< Member to hold the solution matrices
  };

} // namespace solvers
} // namespace sim

#endif // CARLEMAN_SOLVER_HPP
