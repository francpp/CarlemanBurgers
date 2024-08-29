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

  class CarlemanSolver : public Solver
  {
  public:
    CarlemanSolver(
      const params::SimulationParameters          &params,
      const discretization::Discretization        &discretization,
      const initial_conditions::InitialConditions &initialConditions);

    void solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
               Eigen::MatrixXd &F2) override;

    const std::vector<Eigen::MatrixXd> &getUsCN() const;

  private:
    const params::SimulationParameters          &params;
    const discretization::Discretization        &discretization;
    const initial_conditions::InitialConditions &initialConditions;

    Eigen::SparseMatrix<double> prepareCarlemanMatrix(Eigen::MatrixXd &F0,
                                                      Eigen::MatrixXd &F1,
                                                      Eigen::MatrixXd &F2);

    std::vector<Eigen::MatrixXd> us_c_N; // Member to hold the solution matrices
  };

} // namespace solvers
} // namespace sim

#endif // CARLEMAN_SOLVER_HPP
