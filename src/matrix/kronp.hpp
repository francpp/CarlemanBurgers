#ifndef KRONP_HPP
#define KRONP_HPP

#include <Eigen/Sparse>

namespace sim
{
namespace matrix
{

  // Function to calculate the Kronecker product of two matrices using Eigen
  Eigen::SparseMatrix<double> kron(const Eigen::SparseMatrix<double> &A,
                                   const Eigen::SparseMatrix<double> &B);

  // Function to calculate the Kronecker power of a matrix using Eigen
  Eigen::SparseMatrix<double> kronp(const Eigen::SparseMatrix<double> &A,
                                    int                                k);

} // namespace matrix
} // namespace sim

#endif // KRONP_HPP
