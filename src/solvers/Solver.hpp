#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <Eigen/Dense>
#include <vector>

namespace sim::solvers
{
class Solver
{
public:
  virtual ~Solver() = default;
  virtual void solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                     Eigen::MatrixXd &F2) = 0;

  // Common interpolation method
  Eigen::VectorXd interp1(const std::vector<double> &ts,
                          const Eigen::MatrixXd &F0, double t) const;

protected:
  Solver() = default;
};
} // namespace sim::solvers

#endif // SOLVER_HPP
