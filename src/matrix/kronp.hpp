#ifndef KRONP_HPP
#define KRONP_HPP

#include <Eigen/Dense>

namespace sim
{
namespace matrix
{

  // Function to calculate the Kronecker product of two matrices using Eigen
  Eigen::MatrixXd kron(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);

  // Function to calculate the Kronecker power of a matrix using Eigen
  Eigen::MatrixXd kronp(const Eigen::MatrixXd &A, int k);

} // namespace matrix
} // namespace sim

#endif // KRONP_HPP
