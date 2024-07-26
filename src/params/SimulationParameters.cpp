#include "SimulationParameters.hpp"

namespace sim::params
{

std::ostream &
operator<<(std::ostream &out, const SimulationParameters &p)
{
  out << "PARAMETER VALUES:"
      << "\n";
  out << "nx = " << p.nx << "\n";
  out << "nt = " << p.nt << "\n";
  out << "nx_pde = " << p.nx_pde << "\n";
  out << "nt_pde = " << p.nt_pde << "\n";
  out << "Re0 = " << p.Re0 << "\n";
  out << "L0 = " << p.L0 << "\n";
  out << "U0 = " << p.U0 << "\n";
  out << "beta = " << p.beta << "\n";
  out << "f = " << p.f << "\n";
  out << "T = " << p.T << "\n";
  // out << "F0_fun = " << p.F0_fun << "\n";
  out << "N_max = " << p.N_max << "\n";
  out << "ode_deg = " << p.ode_deg << "\n";
  return out;
}

} // namespace sim::params
