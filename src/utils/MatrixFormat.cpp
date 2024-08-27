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

} // namespace matrixUtils
} // namespace sim