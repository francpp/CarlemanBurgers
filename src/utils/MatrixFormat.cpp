#include "MatrixFormat.hpp"

namespace sim
{
namespace matrixUtils
{
  Eigen::MatrixXd
  convertToDenseEigen(const std::vector<std::vector<double>> &vec)
  {
    // Get the dimensions of the input vector
    int rows = vec.size();
    int cols = vec[0].size();

    // Create an Eigen matrix with the same dimensions
    Eigen::MatrixXd eigenMatrix(rows, cols);

    // Copy the data from the 2D vector to the Eigen matrix
    for(int i = 0; i < rows; ++i)
      {
        for(int j = 0; j < cols; ++j)
          {
            eigenMatrix(i, j) = vec[i][j];
          }
      }

    return eigenMatrix;
  }

  void
  assignSparseBlock(Eigen::SparseMatrix<double>       &A,
                    const Eigen::SparseMatrix<double> &Aij, int a0, int b0)
  {
    // Ensure the matrix is in a writable state
    A.reserve(A.nonZeros() + Aij.nonZeros());

    for(int i = 0; i < Aij.outerSize(); ++i)
      {
        for(Eigen::SparseMatrix<double>::InnerIterator it(Aij, i); it; ++it)
          {
            int    row = a0 + it.row();
            int    col = b0 + it.col();
            double value = it.value();

            // Overwrite the value in A with the new value from Aij
            A.coeffRef(row, col) = value;
          }
      }

    // Optionally compress the matrix after modifications
    A.makeCompressed();
  }

} // namespace matrixUtils
} // namespace sim