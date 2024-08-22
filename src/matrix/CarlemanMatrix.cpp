#include "CarlemanMatrix.hpp"
#include "kronp.hpp"
#include <cmath>
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
        dNs[N - 1] = (std::pow(nx, N + 1) - nx) / (nx - 1);
      }

    return dNs;
  }

  Eigen::SparseMatrix<double>
  initializeCarlemanMatrix(int size)
  {
    return Eigen::SparseMatrix<double>(size, size);
  }

  Eigen::SparseMatrix<double>
  assembleCarlemanMatrix(const std::vector<int> &dNs, int N_max, int nx,
                         int ode_deg, const Eigen::SparseMatrix<double> &F0,
                         const Eigen::MatrixXd &F1, const Eigen::MatrixXd &F2)
  {
    int                         total_size = dNs.back();
    Eigen::SparseMatrix<double> A(total_size, total_size);

    // Initialize a full matrix to hold the final Carleman matrix (in dense form
    // first)
    std::vector<std::vector<double>> A_dense(
      total_size, std::vector<double>(total_size, 0.0));

    // Combine F0, F1, F2 into a single matrix-like structure
    std::vector<std::vector<double>> Fs(
      F1.rows(), std::vector<double>(F1.cols() + F2.cols(), 0.0));
    for(int i = 0; i < F1.rows(); ++i)
      {
        for(int j = 0; j < F1.cols(); ++j)
          {
            Fs[i][j] = F1(i, j);
          }
        for(int j = 0; j < F2.cols(); ++j)
          {
            Fs[i][F1.cols() + j] = F2(i, j);
          }
      }

    for(int i = 1; i <= N_max; ++i)
      {
        for(int j = 0; j <= std::min(ode_deg, N_max - i + 1); ++j)
          {
            if(i == 1 && j == 0)
              continue;

            int a0 = (std::pow(nx, i) - nx) / (nx - 1);
            int a1 = a0 + std::pow(nx, i) - 1;
            int b0 = (std::pow(nx, j + i - 1) - nx) / (nx - 1);
            int b1 = b0 + std::pow(nx, j + i - 1) - 1;

            // Initialize Aij to be the correct size, but filled with zeros
            std::vector<std::vector<double>> Aij(
              std::pow(nx, i),
              std::vector<double>(std::pow(nx, j + i - 1), 0.0));

            // Extract the relevant block from Fs
            std::vector<std::vector<double>> Fj(
              nx, std::vector<double>(b1 - b0 + 1, 0.0));
            for(int row = 0; row < nx; ++row)
              {
                for(int col = 0; col < b1 - b0 + 1; ++col)
                  {
                    Fj[row][col] = Fs[row][b0 + col];
                  }
              }

            for(int p = 1; p <= i; ++p)
              {
                auto Ia = kronp(std::vector<std::vector<double>>(
                                  nx, std::vector<double>(nx, 1.0)),
                                p - 1);
                auto Ib = kronp(std::vector<std::vector<double>>(
                                  nx, std::vector<double>(nx, 1.0)),
                                i - p);

                // Perform Kronecker products manually using the provided kron
                // function
                auto kron_product = kron(Ia, Fj);
                for(size_t row = 0; row < Aij.size(); ++row)
                  {
                    for(size_t col = 0; col < Aij[0].size(); ++col)
                      {
                        Aij[row][col] += kron_product[row][col];
                      }
                  }
              }

            // Insert the Aij block into the appropriate place in the full
            // matrix A_dense
            for(int row = 0; row < Aij.size(); ++row)
              {
                for(int col = 0; col < Aij[0].size(); ++col)
                  {
                    A_dense[a0 + row][b0 + col] = Aij[row][col];
                  }
              }
          }
      }

    // Convert A_dense into Eigen::SparseMatrix
    for(int i = 0; i < total_size; ++i)
      {
        for(int j = 0; j < total_size; ++j)
          {
            if(A_dense[i][j] != 0)
              {
                A.insert(i, j) = A_dense[i][j];
              }
          }
      }

    return A;
  }

} // namespace matrix
} // namespace sim
