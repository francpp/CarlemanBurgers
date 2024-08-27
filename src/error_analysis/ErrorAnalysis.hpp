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
class ErrorAnalysis
{
public:
  ErrorAnalysis(const solvers::EulerSolver         &eulerSolver,
                const solvers::ODE45Solver         &odeSolver,
                const solvers::PDESolver           &pdeSolver,
                const solvers::CarlemanSolver      &carlemanSolver,
                const params::SimulationParameters &params);

  void computeErrors();

  // Getters for accessing computed error metrics
  const Eigen::MatrixXd &getEpsCDError() const;
  const Eigen::MatrixXd &getEpsRelCDError() const;
  const Eigen::MatrixXd &getEpsCPDEError() const;
  const Eigen::MatrixXd &getEpsRelCPDEError() const;
  const Eigen::VectorXd &getEpsDPDEError() const;
  const Eigen::VectorXd &getEpsRelDPDEError() const;
  const Eigen::VectorXd &getEpsDEError() const;
  const Eigen::VectorXd &getEpsRelDEError() const;

private:
  // Errors between Carleman solution and ODE solution
  Eigen::MatrixXd eps_c_d_N;     // L2 norm of absolute error
  Eigen::MatrixXd eps_rel_c_d_N; // L-inf norm of relative error

  // Errors between Carleman solution and PDE solution
  Eigen::MatrixXd eps_c_pde_N;     // L2 norm of absolute error
  Eigen::MatrixXd eps_rel_c_pde_N; // L-inf norm of relative error

  // Errors between ODE and PDE solutions
  Eigen::VectorXd eps_d_pde;     // L2 norm of absolute error
  Eigen::VectorXd eps_rel_d_pde; // L-inf norm of relative error

  // Errors between ODE and Euler solutions
  Eigen::VectorXd eps_d_e;     // L2 norm of absolute error
  Eigen::VectorXd eps_rel_d_e; // L-inf norm of relative error

  const solvers::EulerSolver         &eulerSolver;
  const solvers::ODE45Solver         &odeSolver;
  const solvers::PDESolver           &pdeSolver;
  const solvers::CarlemanSolver      &carlemanSolver;
  const params::SimulationParameters &params;
};
} // namespace sim::error_analysis

#endif // ERROR_ANALYSIS_HPP
