#ifndef CARLEMAN_SOLVER_HPP
#define CARLEMAN_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "params/SimulationParameters.hpp"
#include <vector>

namespace sim::solvers
{
class CarlemanSolver : public Solver
{
public:
  std::vector<std::vector<std::vector<double>>> ys_c_N, us_c_N;
  CarlemanSolver(const params::SimulationParameters          &params,
                 const discretization::Discretization        &discretization,
                 const initial_conditions::InitialConditions &ic);
  void solve() override;

private:
  void prepareCarlemanMatrix();
  void solveCarlemanSystem();
};
} // namespace sim::solvers

#endif // CARLEMAN_SOLVER_HPP
