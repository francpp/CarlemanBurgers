#ifndef CARLEMAN_MATRIX_HPP
#define CARLEMAN_MATRIX_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace sim
{
namespace matrix
{

  /**
   * @brief Calculates the block sizes for the Carleman matrix.
   *
   * This function computes the size of each block in the Carleman matrix, which
   * depends on the maximum degree \( N_{\text{max}} \) and the number of
   * spatial points \( nx \).
   *
   * @param N_max The maximum degree of Carleman linearization.
   * @param nx The number of spatial points.
   * @return A vector of integers representing the sizes of each block in the
   * Carleman matrix.
   */
  std::vector<int> calculateBlockSizes(int N_max, int nx);

  /**
   * @brief Assembles the Carleman matrix using the given parameters.
   *
   * The Carleman matrix is constructed using the specified maximum degree \(
   * N_{\text{max}} \), the number of spatial points \( nx \), and the degree of
   * the ODE \( \text{ode\_deg} \). The matrices
   * \( F_0 \), \( F_1 \), and \( F_2 \) are used in the construction of the
   * Carleman matrix.
   *
   * @param N_max The maximum degree of Carleman linearization.
   * @param nx The number of spatial points.
   * @param ode_deg The degree of the ODE.
   * @param F0 The forcing matrix \( F_0 \).
   * @param F1 The linear coefficient matrix \( F_1 \).
   * @param F2 The quadratic interaction matrix \( F_2 \).
   * @return The assembled Carleman matrix as a sparse matrix.
   */
  Eigen::SparseMatrix<double> assembleCarlemanMatrix(int N_max, int nx,
                                                     int ode_deg,
                                                     const Eigen::MatrixXd &F0,
                                                     const Eigen::MatrixXd &F1,
                                                     const Eigen::MatrixXd &F2);

} // namespace matrix
} // namespace sim

#endif // CARLEMAN_MATRIX_HPP
