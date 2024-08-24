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

  double calculateCarlemanConvergenceNumber(const Eigen::MatrixXd     &F0,
                                            const Eigen::MatrixXd     &F1,
                                            const Eigen::MatrixXd     &F2,
                                            const std::vector<double> &u0s,
                                            double dt, int nt, int nx,
                                            int N_max);

} // namespace utils
} // namespace sim

#endif // CARLEMAN_UTILS_HPP
