#include "CarlemanUtils.hpp"

namespace sim
{
namespace utils
{

  // Function to compute the eigenvalues of a matrix (placeholder).
  std::vector<double>
  eig(const std::vector<std::vector<double>> &matrix)
  {
    // Placeholder: return some eigenvalues for the sake of example.
    // TODO
    return {-11.4920, -11.1167, -10.5094, -9.6968, -8.7142, -7.6047,
            -6.4167,  -5.2022,  -4.0142,  -0.1270, -0.5023, -1.1095,
            -1.9222,  -2.9047,  0.0,      0.0};
  }

  // Function to compute the norm of a vector (2-norm).
  double
  norm(const std::vector<double> &vec)
  {
    double sum_of_squares =
      std::accumulate(vec.begin(), vec.end(), 0.0,
                      [](double sum, double val) { return sum + val * val; });
    return std::sqrt(sum_of_squares);
  }

  // Function to compute the norm of a matrix (2-norm, row-wise maximum norm).
  double
  norm(const std::vector<std::vector<double>> &matrix)
  {
    // TODO
    double sum = 0.0;
    for(const auto &row : matrix)
      {
        for(const auto &element : row)
          {
            sum += element * element;
          }
      }
    return std::sqrt(sum);
  }

  double
  calculateCarlemanConvergenceNumber(const std::vector<std::vector<double>> &F0,
                                     const std::vector<std::vector<double>> &F1,
                                     const std::vector<std::vector<double>> &F2,
                                     const std::vector<double> &u0s, double dt,
                                     int nt, int nx, int N_max)
  {
    // Compute eigenvalues of F1
    std::vector<double> lambdas = eig(F1);

    // Remove zero eigenvalues
    lambdas.erase(std::remove_if(lambdas.begin(), lambdas.end(),
                                 [](double lambda) { return lambda == 0; }),
                  lambdas.end());

    // Find the maximum eigenvalue
    double lambda = *std::max_element(lambdas.begin(), lambdas.end());

    std::cout << lambda << std::endl;

    // Compute norms
    double f2 = norm(F2);
    double f1 = norm(F1);
    double f0 = 0.0;
    for(int it = 0; it < nt; ++it)
      {
        f0 = std::max(norm(F0[it]), f0);
      }
    std::cout << f0 << std::endl;
    std::cout << f1 << std::endl;
    std::cout << f2 << std::endl;

    // Calculate R
    double R = (norm(u0s) * f2 + f0 / norm(u0s)) / std::abs(lambda);

    // Calculate r1 and r2
    double discriminant = std::sqrt(lambda * lambda - 4 * f2 * f0);
    double r1 = (std::abs(lambda) - discriminant) / (2 * f2);
    double r2 = (std::abs(lambda) + discriminant) / (2 * f2);

    // Check time step constraint
    if(dt > 1.0 / (N_max * f1))
      {
        throw std::runtime_error("Time step too large");
      }

    // Check perturbation constraint
    if(f0 + f2 > std::abs(lambda))
      {
        std::cout << "Perturbation too large" << std::endl;
      }

    // Returning the calculated R
    return R;
  }

} // namespace utils
} // namespace sim
