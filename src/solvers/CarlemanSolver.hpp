#ifndef CARLEMAN_SOLVER_HPP
#define CARLEMAN_SOLVER_HPP

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

  class CarlemanSolver
  {
  public:
    CarlemanSolver(
      const params::SimulationParameters          &params,
      const discretization::Discretization        &discretization,
      const initial_conditions::InitialConditions &initialConditions);

    void solveCarlemanSystem();

  private:
    const params::SimulationParameters          &params;
    const discretization::Discretization        &discretization;
    const initial_conditions::InitialConditions &initialConditions;

    std::vector<Eigen::SparseMatrix<double>>
                                ys_c_N; // Solution for each truncation level
    Eigen::SparseMatrix<double> A;      // Carleman matrix
    Eigen::VectorXd             b_N;    // Inhomogeneous term
  };

} // namespace solvers
} // namespace sim

#endif // CARLEMAN_SOLVER_HPP
