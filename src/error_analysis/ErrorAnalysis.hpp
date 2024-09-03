#ifndef ERROR_ANALYSIS_HPP
#define ERROR_ANALYSIS_HPP

#include "solvers/CarlemanSolver.hpp"
#include "solvers/EulerSolver.hpp"
#include "solvers/ODE45Solver.hpp"
#include "solvers/PDESolver.hpp"
#include <Eigen/Dense>
#include <vector>

namespace sim::error_analysis
{
/**
 * @class ErrorAnalysis
 * @brief Analyzes the errors between different solvers' solutions.
 *
 * The ErrorAnalysis class computes various error metrics between the solutions
 * obtained from different numerical solvers, such as Carleman, ODE45, PDE, and
 * Euler solvers.
 */
class ErrorAnalysis
{
public:
  /**
   * @brief Constructs an ErrorAnalysis object with references to solvers and
   * parameters.
   * @param eulerSolver Reference to the Euler solver.
   * @param odeSolver Reference to the ODE45 solver.
   * @param pdeSolver Reference to the PDE solver.
   * @param carlemanSolver Reference to the Carleman solver.
   * @param params Reference to the simulation parameters.
   */
  ErrorAnalysis(const solvers::EulerSolver         &eulerSolver,
                const solvers::ODE45Solver         &odeSolver,
                const solvers::PDESolver           &pdeSolver,
                const solvers::CarlemanSolver      &carlemanSolver,
                const params::SimulationParameters &params);

  /**
   * @brief Computes the errors between the solutions of different solvers.
   */
  void computeErrors();

  // Getters for accessing computed error metrics

  /**
   * @brief Retrieves the L2 norm of the absolute error between Carleman and ODE
   * solutions.
   * @return A matrix containing the L2 norm of the absolute error.
   */
  const Eigen::MatrixXd &getEpsCDError() const;

  /**
   * @brief Retrieves the L-inf norm of the relative error between Carleman and
   * ODE solutions.
   * @return A matrix containing the L-inf norm of the relative error.
   */
  const Eigen::MatrixXd &getEpsRelCDError() const;

  /**
   * @brief Retrieves the L2 norm of the absolute error between Carleman and PDE
   * solutions.
   * @return A matrix containing the L2 norm of the absolute error.
   */
  const Eigen::MatrixXd &getEpsCPDEError() const;

  /**
   * @brief Retrieves the L-inf norm of the relative error between Carleman and
   * PDE solutions.
   * @return A matrix containing the L-inf norm of the relative error.
   */
  const Eigen::MatrixXd &getEpsRelCPDEError() const;

  /**
   * @brief Retrieves the L2 norm of the absolute error between ODE and PDE
   * solutions.
   * @return A vector containing the L2 norm of the absolute error.
   */
  const Eigen::VectorXd &getEpsDPDEError() const;

  /**
   * @brief Retrieves the L-inf norm of the relative error between ODE and PDE
   * solutions.
   * @return A vector containing the L-inf norm of the relative error.
   */
  const Eigen::VectorXd &getEpsRelDPDEError() const;

  /**
   * @brief Retrieves the L2 norm of the absolute error between ODE and Euler
   * solutions.
   * @return A vector containing the L2 norm of the absolute error.
   */
  const Eigen::VectorXd &getEpsDEError() const;

  /**
   * @brief Retrieves the L-inf norm of the relative error between ODE and Euler
   * solutions.
   * @return A vector containing the L-inf norm of the relative error.
   */
  const Eigen::VectorXd &getEpsRelDEError() const;

private:
  // Error metrics

  Eigen::MatrixXd eps_c_d_N; ///< L2 norm of absolute error between Carleman and
                             ///< ODE solutions.
  Eigen::MatrixXd eps_rel_c_d_N; ///< L-inf norm of relative error between
                                 ///< Carleman and ODE solutions.
  Eigen::MatrixXd eps_c_pde_N;   ///< L2 norm of absolute error between Carleman
                                 ///< and PDE solutions.
  Eigen::MatrixXd eps_rel_c_pde_N; ///< L-inf norm of relative error between
                                   ///< Carleman and PDE solutions.
  Eigen::VectorXd
    eps_d_pde; ///< L2 norm of absolute error between ODE and PDE solutions.
  Eigen::VectorXd eps_rel_d_pde; ///< L-inf norm of relative error between ODE
                                 ///< and PDE solutions.
  Eigen::VectorXd
    eps_d_e; ///< L2 norm of absolute error between ODE and Euler solutions.
  Eigen::VectorXd eps_rel_d_e; ///< L-inf norm of relative error between ODE and
                               ///< Euler solutions.

  // References to solvers and parameters

  const solvers::EulerSolver         &eulerSolver;
  const solvers::ODE45Solver         &odeSolver;
  const solvers::PDESolver           &pdeSolver;
  const solvers::CarlemanSolver      &carlemanSolver;
  const params::SimulationParameters &params;
};

} // namespace sim::error_analysis

#endif // ERROR_ANALYSIS_HPP
