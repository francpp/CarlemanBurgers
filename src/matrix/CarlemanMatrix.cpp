#include "CarlemanMatrix.hpp"
#include "kronp.hpp"
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace sim
{
namespace matrix
{

  std::vector<int>
  calculateBlockSizes(int N_max, int nx)
  {
    std::vector<int> dNs(N_max);
    for(int N = 1; N <= N_max; ++N)
      {
        int size = static_cast<int>((std::pow(nx, N + 1) - nx) / (nx - 1));
        dNs[N - 1] = size;
      }
    return dNs;
  }

  Eigen::MatrixXd
  assembleCarlemanMatrix(int N_max, int nx, int ode_deg,
                         const Eigen::MatrixXd &F0, const Eigen::MatrixXd &F1,
                         const Eigen::MatrixXd &F2)
  {
    std::vector<int> dNs = calculateBlockSizes(N_max, nx);
    int             total_size = dNs.back();
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(total_size, total_size);

    // Construct Fs by concatenating the first row of F0 as a column vector, F1,
    // and F2
    Eigen::MatrixXd Fs(F1.rows(), 1 + F1.cols() + F2.cols());
    Eigen::MatrixXd firstColumn =
      F0.row(0).transpose(); // Convert the first row of F0 into a column vector
    Fs << firstColumn, F1, F2;
    // TODO
    for(int i = 1; i <= N_max; ++i)
      {
        for (int j = 0; j <= std::min(ode_deg, N_max - i + 1); ++j)
        {
            if (i == 1 && j == 0)
                continue;

            int a0 = 1 + (std::pow(nx, i) - nx) / (nx - 1) -
                     1; // Adjusted for zero-based indexing
            int a1 = a0 + std::pow(nx, i) - 1;
            int b0 = 1 + (std::pow(nx, j + i - 1) - nx) / (nx - 1) -
                     1; // Adjusted for zero-based indexing
            int b1 = b0 + std::pow(nx, j + i - 1) - 1;

            // Initialize Aij using Eigen::MatrixXd
            Eigen::MatrixXd Aij =
              Eigen::MatrixXd::Zero(std::pow(nx, i), std::pow(nx, j + i - 1));
            // Extract the relevant block from Fs
            int f0 = 1 + (std::pow(nx, j) - nx) / (nx - 1) +
                     1; // Adjusted for zero-based indexing
            int f1 = f0 + std::pow(nx, j) -
                     1; // fj and fj+1 difference is the number of rows of Fj
            Eigen::MatrixXd Fj = Fs.block(0, f0 - 1, Fs.rows(), f1 - f0 + 1);
            for (int p = 1; p <= i; ++p)
            {
                // Calculate Kronecker products using Eigen matrices
                Eigen::MatrixXd Ia =
                  kronp(Eigen::MatrixXd::Identity(nx, nx), p - 1);
                Eigen::MatrixXd Ib =
                  kronp(Eigen::MatrixXd::Identity(nx, nx), i - p);

                Eigen::MatrixXd kron_product = kron(Ia, Fj);
                Aij += kron(kron_product,
                            Ib); // Updated based on the MATLAB code's structure
            }

            // Insert the Aij block into the appropriate place in the full
            // matrix A
            A.block(a0, b0, Aij.rows(), Aij.cols()) =
              Aij; // Direct block assignment
        }
      }

    return A;
  }

} // namespace matrix
} // namespace sim
