#include "CarlemanSolver.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include "matrix/kronp.hpp"
#include <iostream>

namespace sim
{
namespace solvers
{

  CarlemanSolver::CarlemanSolver(
    const params::SimulationParameters          &params,
    const discretization::Discretization        &discretization,
    const initial_conditions::InitialConditions &initialConditions)
    : params(params), discretization(discretization),
      initialConditions(initialConditions)
  {}

  void
  CarlemanSolver::solveCarlemanSystem()
  {
    int    nx = params.nx;
    int    nt = params.nt;
    int    N_max = params.N_max;
    double dt = discretization.getTs()[1] - discretization.getTs()[0];
    auto   xs = discretization.getXs();
    auto   ts = discretization.getTs();

    // Create the block sizes
    std::vector<int> dNs = matrix::calculateBlockSizes(N_max, nx);

    // Initialize solution storage
    ys_c_N.resize(N_max);
    for(int N = 0; N < N_max; ++N)
      {
        ys_c_N[N].resize(nt, dNs[N]);
      }

    // Iterate over truncation levels
    for(int N = 1; N <= N_max; ++N)
      {
        Eigen::SparseMatrix<double> A_N = A.block(0, 0, dNs[N - 1], dNs[N - 1]);
        b_N = Eigen::VectorXd::Zero(dNs[N - 1]);

        // Initialize b_N and y0s
        for(int i = 0; i < nx; ++i)
          {
            b_N[i] = initialConditions.getF0()[0][i];
          }

        std::vector<double> y0s;
        for(int i = 1; i <= N; ++i)
          {
            // Wrap u0s in a 2D vector to use with kronp
            std::vector<std::vector<double>> u0s_2D(1,
                                                    initialConditions.getU0s());
            auto kron_product = matrix::kronp(u0s_2D, i);
            y0s.insert(y0s.end(), kron_product[0].begin(),
                       kron_product[0].end());
          }

        std::cout << "Solving Carleman N=" << N << std::endl;

        // Initialize solution matrix ys
        Eigen::MatrixXd ys = Eigen::MatrixXd::Zero(nt, dNs[N - 1]);
        for(int i = 0; i < y0s.size(); ++i)
          {
            ys(0, i) = y0s[i];
          }

        // Time-stepping loop
        for(int k = 0; k < nt - 1; ++k)
          {
            for(int i = 2; i <= N; ++i)
              {
                int a0 = (std::pow(nx, i) - nx) / (nx - 1);
                int a1 = a0 + std::pow(nx, i) - 1;
                int b0 = (std::pow(nx, i - 1) - nx) / (nx - 1);
                int b1 = b0 + std::pow(nx, i - 1) - 1;

                Eigen::SparseMatrix<double> Aij = Eigen::SparseMatrix<double>(
                  std::pow(nx, i), std::pow(nx, i - 1));

                std::vector<std::vector<double>> Fj(
                  nx, std::vector<double>(b1 - b0 + 1, 0.0));
                for(int row = 0; row < nx; ++row)
                  {
                    for(int col = 0; col < b1 - b0 + 1; ++col)
                      {
                        Fj[row][col] = initialConditions.getF0()[k][row];
                      }
                  }

                for(int p = 1; p <= i; ++p)
                  {
                    auto Ia = matrix::kronp(std::vector<std::vector<double>>(
                                              nx, std::vector<double>(nx, 1.0)),
                                            p - 1);
                    auto Ib = matrix::kronp(std::vector<std::vector<double>>(
                                              nx, std::vector<double>(nx, 1.0)),
                                            i - p);

                    auto kron_product = matrix::kron(Ia, Fj);
                    for(int row = 0; row < kron_product.size(); ++row)
                      {
                        for(int col = 0; col < kron_product[0].size(); ++col)
                          {
                            Aij.insert(row, col) = kron_product[row][col];
                          }
                      }
                  }

                A_N.block(a0, b0, a1 - a0 + 1, b1 - b0 + 1) = Aij;
              }

            for(int i = 0; i < nx; ++i)
              {
                b_N[i] = initialConditions.getF0()[k][i];
              }

            ys.row(k + 1) =
              ys.row(k) + dt * (A_N * ys.row(k).transpose() + b_N).transpose();
          }

        std::cout << "Done with N=" << N << std::endl;

        ys_c_N[N - 1] = ys.sparseView();
      }
  }

} // namespace solvers
} // namespace sim
