#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <Eigen/Dense>
#include <vector>

namespace sim::solvers
{
/**
 * @class Solver
 * @brief Abstract base class for solvers in the simulation.
 *
 * The Solver class provides an interface for different types of solvers used
 * in the simulation. It includes a pure virtual function `solve` that must be
 * implemented by derived classes and a common interpolation method `interp1`.
 */
class Solver
{
public:
  /**
   * @brief Virtual destructor for the Solver class.
   */
  virtual ~Solver() = default;

  /**
   * @brief Pure virtual function to solve the system.
   *
   * This function must be implemented by any derived class.
   *
   * @param F0 Forcing matrix F0.
   * @param F1 Linear coefficient matrix F1.
   * @param F2 Quadratic interaction matrix F2.
   */
  virtual void solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                     Eigen::MatrixXd &F2) = 0;

  /**
   * @brief Common interpolation method.
   *
   * Interpolates the value of `F0` at a given time `t` using linear
   * interpolation based on the time vector `ts`.
   *
   * @param ts Vector of time points corresponding to the rows of F0.
   * @param F0 Matrix where each row corresponds to a time point in `ts`.
   * @param t The time at which interpolation is required.
   * @return Interpolated row of F0 as a vector.
   * @throws std::out_of_range if `t` is outside the range of `ts`.
   */
  Eigen::VectorXd interp1(const std::vector<double> &ts,
                          const Eigen::MatrixXd &F0, double t) const;

protected:
  /**
   * @brief Protected default constructor for the Solver class.
   */
  Solver() = default;
};
} // namespace sim::solvers

#endif // SOLVER_HPP
