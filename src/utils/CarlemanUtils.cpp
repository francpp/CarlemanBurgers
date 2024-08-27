#include "CarlemanUtils.hpp"
#include <Eigen/Dense>

namespace sim
{
namespace utils
{

  // Function to compute the eigenvalues of a matrix (placeholder).
  Eigen::VectorXcd
  eig(const Eigen::MatrixXd &matrix)
  {
    Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
    Eigen::VectorXcd                    eigenvalues = solver.eigenvalues();
    return eigenvalues;
  }

  double
  spectralNorm(const Eigen::MatrixXd &matrix)
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    return svd.singularValues()(0); // The largest singular value
  }

  double
  calculateCarlemanConvergenceNumber(const Eigen::MatrixXd     &F0,
                                     const Eigen::MatrixXd     &F1,
                                     const Eigen::MatrixXd     &F2,
                                     const std::vector<double> &u0s, double dt,
                                     int nt, int nx, int N_max)
  {
    // Calculate the eigenvalues of F1
    Eigen::VectorXcd lambdas = eig(F1);

    // Remove zero eigenvalues
    Eigen::VectorXcd nonZeroLambdas;
    for(int i = 0; i < lambdas.size(); ++i)
      {
        if(std::abs(lambdas[i]) > 1e-12)
          { // Small threshold to avoid numerical zero
            nonZeroLambdas.conservativeResize(nonZeroLambdas.size() + 1);
            nonZeroLambdas(nonZeroLambdas.size() - 1) = lambdas[i];
          }
      }

    // Find the maximum eigenvalue
    double lambda = nonZeroLambdas.real().maxCoeff();
    // Calculate norms
    double f2 = spectralNorm(F2);
    double f1 = spectralNorm(F1);
    double f0 = 0.0;

    // Convert F0_vec to Eigen::MatrixXd for the norm calculation

    for(int it = 0; it < nt; ++it)
      {
        double currentNorm = F0.row(it).norm();
        if(currentNorm > f0)
          {
            f0 = currentNorm;
          }
      }

    // Calculate the Carleman convergence number R
    double u0sNorm =
      Eigen::Map<const Eigen::VectorXd>(u0s.data(), u0s.size()).norm();
    double R = (u0sNorm * f2 + f0 / u0sNorm) / std::abs(lambda);

    // Calculate r1 and r2
    double discriminant = lambda * lambda - 4 * f2 * f0;
    double r1 = (std::abs(lambda) - std::sqrt(discriminant)) / (2 * f2);
    double r2 = (std::abs(lambda) + std::sqrt(discriminant)) / (2 * f2);

    // Check conditions and print warnings or errors
    if(dt > 1 / (N_max * f1))
      {
        throw std::runtime_error("Time step too large");
      }

    if(f0 + f2 > std::abs(lambda))
      {
        std::cerr << "Perturbation too large" << std::endl;
      }

    return R;
  }

} // namespace utils
} // namespace sim
