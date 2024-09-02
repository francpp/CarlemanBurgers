#include "InitialConditions.hpp"
#include "utils/muparser_fun.hpp"
#include <algorithm>
#include <cmath>

namespace sim::initial_conditions
{

InitialConditions::InitialConditions(
  const sim::params::SimulationParameters   &params,
  const sim::discretization::Discretization &discretization)
  : params(params), discretization(discretization)
{}

void
InitialConditions::computeInitialConditions()
{
  const auto &xs = discretization.getXs();
  u0s.resize(xs.size());

  // Compute initial condition u0(x) = -U0 * sin(2 * pi * f * x / L0)
  double x_var = 0.0; // placeholder for x
  double U0 = params.U0;
  double L0 = params.L0;
  double f = params.f;

  std::string                     u0_expr = "-U0*sin(2*pi*f*x/L0)";
  std::map<std::string, double *> u0_vars = {
    {"x", &x_var}, {"U0", &U0}, {"L0", &L0}, {"f", &f}};
  utils::MuparserFun u0_fun(u0_expr, u0_vars);

  std::transform(xs.begin(), xs.end(), u0s.begin(), [&](double x) {
    std::map<std::string, double> u0_values = {{"x", x}};
    return u0_fun(u0_values);
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
  double      nu = params.nu;

  // Initialize F0
  F0.resize(nt, std::vector<double>(nx, 0));
  double t_var = 0.0; // placeholder for t
  double x_var = 0.0; // placeholder for x
  double U0 = params.U0;
  double L0 = params.L0;

  std::string F0_expr = "U0*exp(-(x-L0/4)^2/(2*(L0/32)^2))*cos(2*pi*t)";
  std::map<std::string, double *> F0_vars = {
    {"t", &t_var}, {"x", &x_var}, {"U0", &U0}, {"L0", &L0}};
  utils::MuparserFun F0_fun(F0_expr, F0_vars);

  for(int it = 0; it < nt; ++it)
    {
      std::transform(xs.begin(), xs.end(), F0[it].begin(), [&](double x) {
        std::map<std::string, double> F0_values = {{"t", ts[it]}, {"x", x}};

        return F0_fun(F0_values);
      });
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
  for(int i = 0; i < nx; ++i) // Loop over rows
    {
      size_t index1 = (i - 1) * nx + i - 1;
      size_t index2 = (i + 1) * nx + i + 1;
      if((index1 < nx * nx) && (index2 < nx * nx))
        {
          F2[i][index1] = 1 / (4 * dx);
          F2[i][index2] = -1 / (4 * dx);
        }
    }
  // Enforce Dirichlet boundaries within the domain
  std::fill(F2[0].begin(), F2[0].end(), 0);
  std::fill(F2[nx - 1].begin(), F2[nx - 1].end(), 0);
}

// Overload for printing Discretization parameters, including the new ones
std::ostream &
operator<<(std::ostream &out, const InitialConditions &d)
{
  out << "Initial Condition PARAMETERS:" << std::endl;
  out << "-------------------------" << std::endl;
  // print F0, F1, F2, u0s
  out << "F0: ";
  for(int i = 0; i < std::min(25, static_cast<int>(d.getF0().size())); ++i)
    {
      for(int j = 0; j < std::min(25, static_cast<int>(d.getF0()[i].size()));
          ++j)
        {
          if(std::abs(d.getF0()[i][j]) > 0.0000001)
            out << d.getF0()[i][j] << " ";
          else
            out << "0 ";
        }
      out << std::endl;
    }

  out << "F1: ";
  for(int i = 0; i < std::min(25, static_cast<int>(d.getF1().size())); ++i)
    {
      for(int j = 0; j < std::min(25, static_cast<int>(d.getF1()[i].size()));
          ++j)
        {
          out << d.getF1()[i][j] << " ";
        }
      out << std::endl;
    }

  out << "F2: ";
  for(int i = 0; i < std::min(25, static_cast<int>(d.getF2().size())); ++i)
    {
      for(int j = 0; j < std::min(25, static_cast<int>(d.getF2()[i].size()));
          ++j)
        {
          out << d.getF2()[i][j] << " ";
        }
      out << std::endl;
    }

  out << "u0s: ";
  for(int i = 0; i < std::min(10, static_cast<int>(d.getU0s().size())); ++i)
    out << d.getU0s()[i] << " ";

  out << std::endl;

  return out;
}

} // namespace sim::initial_conditions
