#include "ErrorAnalysis.hpp"

namespace sim::error_analysis
{

ErrorAnalysis::ErrorAnalysis(const sim::solvers::EulerSolver    &eulerSolver,
                             const sim::solvers::ODE45Solver    &odeSolver,
                             const sim::solvers::PDESolver      &pdeSolver,
                             const sim::solvers::CarlemanSolver &carlemanSolver)
{
  // Constructor implementation
}

void
ErrorAnalysis::computeErrors()
{
  // Error computation implementation
}

} // namespace sim::error_analysis
