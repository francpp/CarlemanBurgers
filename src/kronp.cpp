#include <Eigen/Eigen>

Eigen::MatrixXd
kronp(const Eigen::MatrixXd &A, int k)
{
  Eigen::MatrixXd B = Eigen::MatrixXd::Ones(1, 1);
  for(int i = 0; i < k; i++)
    {
      B = Eigen::kroneckerProduct(B, A);
    }
  return B;
}