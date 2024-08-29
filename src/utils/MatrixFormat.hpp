#ifndef MATRIX_FORMAT_HPP
#define MATRIX_FORMAT_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sim
{
namespace matrixUtils
{
  Eigen::MatrixXd convertToDenseEigen(const std::vector<std::vector<double>> &);
  void            assignSparseBlock(Eigen::SparseMatrix<double> &,
                                    const Eigen::SparseMatrix<double> &, int, int);

} // namespace matrixUtils
} // namespace sim

#endif // MATRIX_FORMAT_HPP
