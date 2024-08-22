#ifndef CARLEMAN_UTILS_HPP
#define CARLEMAN_UTILS_HPP

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace sim
{
namespace utils
{

  double
  calculateCarlemanConvergenceNumber(const std::vector<std::vector<double>> &F0,
                                     const std::vector<std::vector<double>> &F1,
                                     const std::vector<std::vector<double>> &F2,
                                     const std::vector<double> &u0s, double dt,
                                     int nt, int nx, int N_max);

} // namespace utils
} // namespace sim

#endif // CARLEMAN_UTILS_HPP
