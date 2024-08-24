#include "kronp.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace sim
{
namespace matrix
{

  // Function to calculate the Kronecker product of two matrices using Eigen
  Eigen::MatrixXd
  kron(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
  {
    const int aRows = A.rows();
    const int aCols = A.cols();
    const int bRows = B.rows();
    const int bCols = B.cols();

    Eigen::MatrixXd C(aRows * bRows, aCols * bCols);

    for(int i = 0; i < aRows; ++i)
      {
        for(int j = 0; j < aCols; ++j)
          {
            C.block(i * bRows, j * bCols, bRows, bCols) = A(i, j) * B;
          }
      }

    return C;
  }

  // Function to calculate the Kronecker power of a matrix using Eigen
  Eigen::MatrixXd
  kronp(const Eigen::MatrixXd &A, int k)
  {
    Eigen::MatrixXd B(1, 1);
    B(0, 0) = 1.0;

    for(int i = 0; i < k; ++i)
      {
        B = kron(B, A);
      }

    return B;
  }

} // namespace matrix
} // namespace sim
