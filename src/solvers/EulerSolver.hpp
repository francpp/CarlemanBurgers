#ifndef EULER_SOLVER_HPP
#define EULER_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include <vector>

namespace sim::solvers
{
class EulerSolver : public Solver
{
public:
  std::vector<std::vector<double>> us_e;
  EulerSolver(const params::SimulationParameters          &params,
              const discretization::Discretization        &discretization,
              const initial_conditions::InitialConditions &ic);
  void solve() override;
};
} // namespace sim::solvers

#endif // EULER_SOLVER_HPP
