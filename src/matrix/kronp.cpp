#include "kronp.hpp"
#include "utils/MatrixFormat.hpp"
#include <Eigen/Sparse>
#include <iostream>
namespace sim
{
namespace matrix
{

  // Function to calculate the Kronecker product of two matrices using Eigen
  // (Sparse version)
  Eigen::SparseMatrix<double>
  kron(const Eigen::SparseMatrix<double> &A,
       const Eigen::SparseMatrix<double> &B)
  {
    const int aRows = A.rows();
    const int aCols = A.cols();
    const int bRows = B.rows();
    const int bCols = B.cols();

    Eigen::SparseMatrix<double> C(aRows * bRows, aCols * bCols);

    for(int k = 0; k < A.outerSize(); ++k)
      {
        for(Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
          {
            matrixUtils::assignSparseBlock(C, it.value() * B, it.row() * bRows,
                                           it.col() * bCols);
          }
      }

    return C;
  }

  // Function to calculate the Kronecker power of a matrix using Eigen (Sparse
  // version)
  Eigen::SparseMatrix<double>
  kronp(const Eigen::SparseMatrix<double> &A, int k)
  {
    Eigen::SparseMatrix<double> B(1, 1);
    B.insert(0, 0) = 1.0;

    for(int i = 0; i < k; ++i)
      {
        B = kron(B, A);
      }

    return B;
  }

} // namespace matrix
} // namespace sim
