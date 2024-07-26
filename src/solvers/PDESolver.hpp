#ifndef PDE_SOLVER_HPP
#define PDE_SOLVER_HPP

#include "Solver.hpp"
#include "discretization/Discretization.hpp"
#include "initial_conditions/InitialConditions.hpp"
#include "matrix/MatrixOperations.hpp"
#include "params/SimulationParameters.hpp"
#include <vector>

namespace sim::solvers
{
class PDESolver : public Solver
{
public:
  std::vector<std::vector<double>> us_pde;
  PDESolver(const params::SimulationParameters          &params,
            const discretization::Discretization        &discretization,
            const initial_conditions::InitialConditions &ic,
            const matrix::MatrixOperations              &mo);
  void solve() override;
};
} // namespace sim::solvers

#endif // PDE_SOLVER_HPP
