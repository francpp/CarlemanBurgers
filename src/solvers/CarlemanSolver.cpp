#include "CarlemanSolver.hpp"
#include "matrix/CarlemanMatrix.hpp"
#include "matrix/kronp.hpp"
#include <iostream>
#include <vector>

namespace sim
{
namespace solvers
{

  CarlemanSolver::CarlemanSolver(
    const params::SimulationParameters          &params,
    const discretization::Discretization        &discretization,
    const initial_conditions::InitialConditions &initialConditions)
    : params(params), discretization(discretization),
      initialConditions(initialConditions), us_c_N(params.N_max)
  {}

  Eigen::MatrixXd
  CarlemanSolver::prepareCarlemanMatrix(Eigen::MatrixXd &F0,
                                        Eigen::MatrixXd &F1,
                                        Eigen::MatrixXd &F2)
  {
    return matrix::assembleCarlemanMatrix(params.N_max, params.nx,
                                          params.ode_deg, F0, F1, F2);
  }

  void
  CarlemanSolver::solve(Eigen::MatrixXd &F0, Eigen::MatrixXd &F1,
                        Eigen::MatrixXd &F2)
  {
    Eigen::MatrixXd carleman_matrix = prepareCarlemanMatrix(F0, F1, F2);
    int             nx = params.nx;
    int             nt = params.nt;
    int             N_max = params.N_max;
    double          dt = discretization.getTs()[1] - discretization.getTs()[0];

    Eigen::MatrixXd xs = Eigen::Map<const Eigen::MatrixXd>(
      discretization.getXs().data(), discretization.getXs().size(), 1);
    Eigen::MatrixXd ts = Eigen::Map<const Eigen::MatrixXd>(
      discretization.getTs().data(), discretization.getTs().size(), 1);
    Eigen::MatrixXd u0s = Eigen::Map<const Eigen::MatrixXd>(
      initialConditions.getU0s().data(), initialConditions.getU0s().size(), 1);
    std::vector<int> dNs = matrix::calculateBlockSizes(N_max, nx);

    double U0 = params.U0;
    double L0 = params.L0;

    std::function<Eigen::MatrixXd(double, const Eigen::MatrixXd &)> F0_fun =
      [U0, L0](double t, const Eigen::MatrixXd &xs) {
        Eigen::MatrixXd gaussian =
          (-((xs.array() - L0 / 4).square()) / (2 * std::pow(L0 / 32, 2)))
            .exp();
        double          cos_component = std::cos(2 * M_PI * t);
        Eigen::MatrixXd result = U0 * gaussian * cos_component;
        return result;
      };

    for(int N = 1; N <= N_max; ++N)
      {
        int             dN = dNs[N - 1];
        Eigen::MatrixXd A_N = carleman_matrix.block(0, 0, dN, dN);
        Eigen::MatrixXd b_N = Eigen::MatrixXd::Zero(dN, 1);

        b_N.block(0, 0, nx, 1) = F0_fun(ts(0), xs);
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

        std::vector<Eigen::MatrixXd> ys(nt);
        ys[0] = y0s;

        for(int k = 0; k < nt - 1; ++k)
          {
            double current_t = ts(k);

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
                    Aij += matrix::kron(matrix::kron(Ia, Fj), Ib).sparseView();
                  }

                A_N.block(a0, b0, a1 - a0 + 1, b1 - b0 + 1) = Aij;
              }

            b_N.block(0, 0, nx, 1) = F0_fun(current_t, xs);
            ys[k + 1] = ys[k] + dt * (A_N * ys[k] + b_N);
          }

        us_c_N[N - 1] = Eigen::MatrixXd(nt, nx);
        for(int k = 0; k < nt; ++k)
          {
            us_c_N[N - 1].row(k) = ys[k].block(0, 0, nx, 1).transpose();
          }
      }
  }

  const std::vector<Eigen::MatrixXd> &
  CarlemanSolver::getUsCN() const
  {
    return us_c_N;
  }

} // namespace solvers
} // namespace sim
