#include "InitialConditions.hpp"
#include <algorithm>
#include <cmath>

namespace sim::initial_conditions
{

InitialConditions::InitialConditions(
  const sim::params::SimulationParameters   &params,
  const sim::discretization::Discretization &discretization)
  : params(params), discretization(discretization)
{
  // Define the source function F0_fun
  F0_fun = [&](double t, double x) {
    return params.U0 *
           std::exp(-std::pow(x - params.L0 / 4, 2) /
                    (2 * std::pow(params.L0 / 32, 2))) *
           std::cos(2 * M_PI * t);
  };
}

void
InitialConditions::computeInitialConditions()
{
  const auto &xs = discretization.getXs();
  u0s.resize(xs.size());

  // Compute initial condition u0(x) = -U0 * sin(2 * pi * f * x / L0)
  std::transform(xs.begin(), xs.end(), u0s.begin(), [&](double x) {
    return -params.U0 * std::sin(2 * M_PI * params.f * x / params.L0);
  });
}

void
InitialConditions::computeForcingBoundaryConditions()
{
  const auto &xs = discretization.getXs();
  const auto &ts = discretization.getTs();
  int         nx = xs.size();
  int         nt = ts.size();
  double      dx = xs[1] - xs[0];
  double      nu = params.U0 * params.L0 / params.Re0; // Viscosity

  // Initialize F0
  F0.resize(nt, std::vector<double>(nx, 0));
  for(int it = 0; it < nt; ++it)
    {
      std::transform(xs.begin(), xs.end(), F0[it].begin(),
                     [&](double x) { return F0_fun(ts[it], x); });
    }

  // Initialize F1
  F1.resize(nx, std::vector<double>(nx, 0));
  for(int i = 1; i < nx - 1; ++i)
    {
      F1[i][i] = -2 * nu / (dx * dx) - params.beta;
      F1[i][i - 1] = nu / (dx * dx);
      F1[i][i + 1] = nu / (dx * dx);
    }
  // Enforce Dirichlet boundaries within the domain
  F1[0] = std::vector<double>(nx, 0);
  F1[nx - 1] = std::vector<double>(nx, 0);

  // Initialize F2
  F2.resize(nx, std::vector<double>(nx * nx, 0));
  for(int i = 1; i < nx - 1; ++i)
    {
      F2[i][i * nx + i - 1] = -1 / (4 * dx);
      F2[i][i * nx + i + 1] = 1 / (4 * dx);
    }
  // Enforce Dirichlet boundaries within the domain
  std::fill(F2[0].begin(), F2[0].end(), 0);
  std::fill(F2[nx - 1].begin(), F2[nx - 1].end(), 0);
}

} // namespace sim::initial_conditions
