#ifndef CARLEMAN_UTILS_HPP
#define CARLEMAN_UTILS_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace sim
{
namespace utils
{
  /**
   * @brief Calculates the Carleman convergence number for a given system.
   *
   * The Carleman convergence number is used to assess the stability and
   * convergence of the Carleman linearization approach applied to a nonlinear
   * system of equations.
   *
   * @param F0 The matrix representing the forcing function.
   * @param F1 The matrix containing the linear coefficients of the system.
   * @param F2 The matrix representing the quadratic interactions in the system.
   * @param u0s The initial condition vector.
   * @param dt The time step size.
   * @param nt The number of time steps.
   * @param nx The number of spatial points.
   * @param N_max The maximum degree of Carleman linearization.
   * @return The Carleman convergence number \(R\).
   * @throws std::runtime_error if the time step is too large.
   */
  double calculateCarlemanConvergenceNumber(const Eigen::MatrixXd     &F0,
                                            const Eigen::MatrixXd     &F1,
                                            const Eigen::MatrixXd     &F2,
                                            const std::vector<double> &u0s,
                                            double dt, int nt, int nx,
                                            int N_max);

} // namespace utils
} // namespace sim

#endif // CARLEMAN_UTILS_HPP
