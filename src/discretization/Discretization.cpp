#include "Discretization.hpp"
#include <cmath>

namespace sim::discretization
{

Discretization::Discretization(const sim::params::SimulationParameters &params)
  : params(params)
{}

void
Discretization::createDiscretization()
{
  // Destructure the necessary parameters from getParams() method
  const auto &[nx, nt, nx_pde, nt_pde, Re0, L0, U0, beta, f, T, N_max,
               ode_deg] = getParams();

  double x0 = -L0 / 2;
  double x1 = L0 / 2;
  double t0 = 0;
  double t1 = T;

  dx = (x1 - x0) / (nx - 1);
  dt = (t1 - t0) / (nt - 1);

  xs.resize(nx);
  ts.resize(nt);

  for(int i = 0; i < nx; ++i)
    xs[i] = x0 + i * dx;
  for(int i = 0; i < nt; ++i)
    ts[i] = t0 + i * dt;

  dx_pde = (x1 - x0) / (nx_pde - 1);
  dt_pde = (t1 - t0) / (nt_pde - 1);

  ts_pde.resize(nt_pde);
  xs_pde.resize(nx_pde);

  for(int i = 0; i < nx_pde; ++i)
    xs_pde[i] = x0 + i * dx_pde;
  for(int i = 0; i < nt_pde; ++i)
    ts_pde[i] = t0 + i * dt_pde;
}

// Overload for printing Discretization parameters, including the new ones
std::ostream &
operator<<(std::ostream &out, const Discretization &d)
{
  out << "DISCRETIZATION PARAMETERS:" << std::endl;
  out << "-------------------------" << std::endl;

  out << "nx: " << d.getXs().size() << std::endl;
  out << "nt: " << d.getTs().size() << std::endl;

  out << "xs: ";
  for(int i = 0; i < std::min(10, static_cast<int>(d.getXs().size())); ++i)
    out << d.getXs()[i] << " ";
  out << std::endl;

  out << "ts: ";
  for(int i = 0; i < std::min(10, static_cast<int>(d.getTs().size())); ++i)
    out << d.getTs()[i] << " ";
  out << std::endl;

  std::cout << "\n";

  out << "nx_pde: " << d.getXsPde().size() << std::endl;
  out << "nt_pde: " << d.getTsPde().size() << std::endl;

  out << "xs_pde: ";
  for(int i = 0; i < std::min(10, static_cast<int>(d.getXsPde().size())); ++i)
    out << d.getXsPde()[i] << " ";
  out << std::endl;

  out << "ts_pde: ";
  for(int i = 0; i < std::min(10, static_cast<int>(d.getTsPde().size())); ++i)
    out << d.getTsPde()[i] << " ";
  out << std::endl;

  return out;
}

} // namespace sim::discretization
