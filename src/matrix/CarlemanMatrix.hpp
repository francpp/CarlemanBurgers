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

  // Function to initialize the Carleman matrix as a sparse matrix
  Eigen::SparseMatrix<double> initializeCarlemanMatrix(int size);

  // Function to assemble the Carleman matrix using the given parameters
  Eigen::SparseMatrix<double>
  assembleCarlemanMatrix(const std::vector<int> &dNs, int N_max, int nx,
                         int ode_deg, const Eigen::SparseMatrix<double> &F0,
                         const Eigen::MatrixXd &F1, const Eigen::MatrixXd &F2);

} // namespace matrix
} // namespace sim

#endif // CARLEMAN_MATRIX_HPP
