#ifndef MATRIX_FORMAT_HPP
#define MATRIX_FORMAT_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace sim
{
namespace matrixUtils
{
  /**
   * @brief Converts a 2D vector of doubles to an Eigen dense matrix.
   *
   * This function takes a standard 2D vector of doubles and converts it into
   * an Eigen `MatrixXd` for further numerical operations.
   *
   * @param vec A 2D vector of doubles representing the matrix.
   * @return An Eigen `MatrixXd` containing the same data as the input vector.
   */
  Eigen::MatrixXd
  convertToDenseEigen(const std::vector<std::vector<double>> &vec);

  /**
   * @brief Assigns a sparse block matrix into a larger sparse matrix at a
   * specified position.
   *
   * This function inserts a smaller sparse matrix `Aij` into a larger sparse
   * matrix `A` starting at the position (`a0`, `b0`).
   *
   * @param A The larger sparse matrix where the block will be assigned.
   * @param Aij The smaller sparse matrix block to assign into `A`.
   * @param a0 The starting row index in `A` where `Aij` will be placed.
   * @param b0 The starting column index in `A` where `Aij` will be placed.
   */
  void assignSparseBlock(Eigen::SparseMatrix<double>       &A,
                         const Eigen::SparseMatrix<double> &Aij, int a0,
                         int b0);

} // namespace matrixUtils
} // namespace sim

#endif // MATRIX_FORMAT_HPP
