#ifndef MATRIX_FORMAT_HPP
#define MATRIX_FORMAT_HPP

#include <Eigen/Dense>

namespace sim
{
namespace matrixUtils
{
  Eigen::MatrixXd convertToDenseEigen(const std::vector<std::vector<double>> &);

} // namespace matrixUtils
} // namespace sim

#endif // MATRIX_FORMAT_HPP
