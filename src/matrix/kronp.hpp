#ifndef KRONP_HPP
#define KRONP_HPP

#include <Eigen/Sparse>

namespace sim
{
namespace matrix
{

  /**
   * @brief Computes the Kronecker product of two sparse matrices.
   *
   * The Kronecker product of matrices \( A \) and \( B \) is a block matrix
   * where each element \( a_{ij} \) of matrix \( A \) is multiplied by the
   * matrix \( B \).
   *
   * @param A The first sparse matrix.
   * @param B The second sparse matrix.
   * @return A sparse matrix representing the Kronecker product of \( A \) and
   * \( B \).
   */
  Eigen::SparseMatrix<double> kron(const Eigen::SparseMatrix<double> &A,
                                   const Eigen::SparseMatrix<double> &B);

  /**
   * @brief Computes the Kronecker power of a sparse matrix.
   *
   * The Kronecker power of a matrix \( A \) with exponent \( k \) is defined as
   * the Kronecker product of \( A \) with itself \( k \) times.
   *
   * @param A The sparse matrix to be raised to the power \( k \) using
   * Kronecker product.
   * @param k The exponent indicating how many times the Kronecker product
   * should be applied.
   * @return A sparse matrix representing the Kronecker power of \( A \).
   */
  Eigen::SparseMatrix<double> kronp(const Eigen::SparseMatrix<double> &A,
                                    int                                k);

} // namespace matrix
} // namespace sim

#endif // KRONP_HPP
