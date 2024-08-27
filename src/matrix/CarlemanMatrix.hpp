#ifndef CARLEMAN_MATRIX_HPP
#define CARLEMAN_MATRIX_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace sim
{
namespace matrix
{
  // Function to calculate matrix block sizes for the Carleman matrix
  std::vector<int> calculateBlockSizes(int N_max, int nx);

  // Function to assemble the Carleman matrix using the given parameters
  Eigen::MatrixXd assembleCarlemanMatrix(int N_max, int nx, int ode_deg,
                                         const Eigen::MatrixXd &F0,
                                         const Eigen::MatrixXd &F1,
                                         const Eigen::MatrixXd &F2);

} // namespace matrix
} // namespace sim

#endif // CARLEMAN_MATRIX_HPP
