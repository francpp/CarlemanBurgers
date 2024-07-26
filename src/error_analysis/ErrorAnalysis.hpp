#ifndef ERROR_ANALYSIS_HPP
#define ERROR_ANALYSIS_HPP

#include "solvers/CarlemanSolver.hpp"
#include "solvers/EulerSolver.hpp"
#include "solvers/ODE45Solver.hpp"
#include "solvers/PDESolver.hpp"
#include <vector>

namespace sim::error_analysis
{
class ErrorAnalysis
{
public:
  std::vector<std::vector<double>> eps_c_d_N, eps_rel_c_d_N, eps_c_pde_N,
    eps_rel_c_pde_N, eps_d_pde, eps_rel_d_pde, eps_d_e, eps_rel_d_e;

  ErrorAnalysis(const solvers::EulerSolver    &eulerSolver,
                const solvers::ODE45Solver    &odeSolver,
                const solvers::PDESolver      &pdeSolver,
                const solvers::CarlemanSolver &carlemanSolver);
  void computeErrors();
};
} // namespace sim::error_analysis

#endif // ERROR_ANALYSIS_HPP
