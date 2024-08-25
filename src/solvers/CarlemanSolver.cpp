#include "CarlemanSolver.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include "matrix/kronp.hpp"
#include <fstream>
#include <iostream>
#include <vector>

namespace sim
{
namespace solvers
{
  // need to pass also the matrice A and dNs
  CarlemanSolver::CarlemanSolver(
    const params::SimulationParameters          &params,
    const discretization::Discretization        &discretization,
    const initial_conditions::InitialConditions &initialConditions)
    : params(params), discretization(discretization),
      initialConditions(initialConditions)
  {}

  void
  CarlemanSolver::solveCarlemanSystem(Eigen::MatrixXd &A)
  {
    int             nx = params.nx;
    int             nt = params.nt;
    int             N_max = params.N_max;
    double dt = discretization.getTs()[1] - discretization.getTs()[0];
    Eigen::MatrixXd xs = Eigen::Map<const Eigen::MatrixXd>(
      discretization.getXs().data(), discretization.getXs().size(), 1);
    Eigen::MatrixXd ts = Eigen::Map<const Eigen::MatrixXd>(
      discretization.getTs().data(), discretization.getTs().size(), 1);
    Eigen::MatrixXd u0s = Eigen::Map<const Eigen::MatrixXd>(
      initialConditions.getU0s().data(), initialConditions.getU0s().size(), 1);
    std::vector<int> dNs = matrix::calculateBlockSizes(N_max, nx);

    // Use U0 and L0 from params
    double U0 = params.U0;
    double L0 = params.L0;

    // Define F0_fun as a lambda function using U0 and L0 from params
    std::function<Eigen::MatrixXd(double, const Eigen::MatrixXd &)> F0_fun =
      [U0, L0](double t, const Eigen::MatrixXd &xs) {
        // Compute the Gaussian component
        Eigen::MatrixXd gaussian =
          (-((xs.array() - L0 / 4).square()) / (2 * std::pow(L0 / 32, 2)))
            .exp();
        // Compute the cosine component
        double cos_component = std::cos(2 * M_PI * t);

        Eigen::MatrixXd result = U0 * gaussian * cos_component;

        // Return the product of U0, the Gaussian, and the cosine component
        return result;
      };
    // Initialize ys_c_N as a vector of vectors of Eigen::MatrixXd to simulate a
    // 3D array
    std::vector<std::vector<Eigen::MatrixXd>> ys_c_N(
      N_max, std::vector<Eigen::MatrixXd>(nt));
    for(int N = 1; N <= N_max; ++N)
      {
        int dN = dNs[N - 1]; // Adjusted for 0-based indexing

        Eigen::MatrixXd A_N = A.block(0, 0, dN, dN);
        Eigen::MatrixXd b_N = Eigen::MatrixXd::Zero(dN, 1);

        // Initialize b_N for the first time step
        b_N.block(0, 0, nx, 1) = F0_fun(ts(0), xs);
        // Construct y0s
        Eigen::MatrixXd y0s = Eigen::MatrixXd::Zero(nt, dN);

        for(int i = 1; i <= N; ++i)
          {
            Eigen::MatrixXd kron_result = matrix::kronp(u0s, i);
            if(i == 1)
              {
                y0s = kron_result;
              }
            else
              {
                Eigen::MatrixXd temp(y0s.rows() + kron_result.rows(), 1);
                temp << y0s, kron_result;
                y0s = temp;
              }
          }

        std::cout << "Solving Carleman N=" << N << std::endl;

        // Initialize ys for all time steps
        std::vector<Eigen::MatrixXd> ys(nt);
        ys[0] = y0s;

        for(int k = 0; k < nt - 1; ++k)
          {
            double current_t = ts(k);

            // Rebuild the inhomogeneous part of the Carleman matrix per time
            // step
            for(int i = 2; i <= N; ++i)
              {
                int a0 =
                  (nx == 1) ? 0 : (1 + (std::pow(nx, i) - nx) / (nx - 1)) - 1;
                int a1 = a0 + std::pow(nx, i) - 1;
                int b0 = (nx == 1)
                           ? 0
                           : (1 + (std::pow(nx, i - 1) - nx) / (nx - 1)) - 1;
                int b1 = b0 + std::pow(nx, i - 1) - 1;

                Eigen::SparseMatrix<double> Aij(std::pow(nx, i),
                                                std::pow(nx, i - 1));
                Eigen::MatrixXd             Fj = F0_fun(current_t, xs);
                for(int p = 1; p <= i; ++p)
                  {
                    Eigen::MatrixXd Ia =
                      matrix::kronp(Eigen::MatrixXd::Identity(nx, nx), p - 1);
                    Eigen::MatrixXd Ib =
                      matrix::kronp(Eigen::MatrixXd::Identity(nx, nx), i - p);
                    Eigen::MatrixXd kron_result =
                      matrix::kron(matrix::kron(Ia, Fj), Ib);
                    Aij += kron_result.sparseView();
                  }

                A_N.block(a0, b0, a1 - a0 + 1, b1 - b0 + 1) = Aij;
              }

            b_N.block(0, 0, nx, 1) = F0_fun(current_t, xs);
            ys[k + 1] = ys[k] + dt * (A_N * ys[k] + b_N);
          }

        // Store results in ys_c_N

        for(int k = 0; k < nt; ++k)
          {
            ys_c_N[N - 1][k] = ys[k].block(0, 0, dN, 1);
          }
      }

    // Extract us_c_N

    std::vector<std::vector<Eigen::MatrixXd>> us_c_N(
      N_max, std::vector<Eigen::MatrixXd>(nt));
    for(int N = 1; N <= N_max; ++N)
      {
        for(int k = 0; k < nt; ++k)
          {
            us_c_N[N - 1][k] = ys_c_N[N - 1][k].block(0, 0, nx, 1);
          }
      }

    std::cout << "Carleman system solved" << std::endl;
    // cout size and first values
    std::cout << "us_c_N size: " << us_c_N.size() << std::endl;
    std::cout << "us_c_N[0] size: " << us_c_N[0].size() << std::endl;
    std::cout << "us_c_N[0][0] size: " << us_c_N[0][0].size() << std::endl;
    std::cout << "us_c_N[0][0] first value: " << us_c_N[2][5] << std::endl;

    // Open the file to write the data
    std::ofstream outfile("output.txt");
    if(!outfile.is_open())
      {
        std::cerr << "Error opening file!" << std::endl;
      }

    // Iterate over the vector and write the data
    for(int i = 0; i < N_max; ++i)
      {
        for(int j = 0; j < nt; ++j)
          {
            outfile << "Matrix (" << i << ", " << j << "):\n";
            outfile << us_c_N[i][j] << "\n\n";
          }
      }
  }

} // namespace solvers
} // namespace sim
