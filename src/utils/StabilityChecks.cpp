#include "StabilityChecks.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace sim
{
namespace stability
{

  void
  checkCFLConditions(double U0, double dt, double dx, double nu, double dt_ode,
                     double dx_pde, double dt_pde)
  {
    double C1_e = U0 * dt / dx;
    double C2_e = 2 * nu * dt / (dx * dx);
    double C1_ode = U0 * dt_ode / dx;
    double C2_ode = 2 * nu * dt_ode / (dx * dx);
    double C1_pde = U0 * dt_pde / dx_pde;
    double C2_pde = 2 * nu * dt_pde / (dx_pde * dx_pde);

    std::ostringstream error_message;

    if(C1_e > 1)
      {
        error_message << "C1_e = " << C1_e << " exceeds stability limit.\n";
      }
    if(C2_e > 1)
      {
        error_message << "C2_e = " << C2_e << " exceeds stability limit.\n";
      }
    if(C1_ode > 1)
      {
        error_message << "C1_ode = " << C1_ode << " exceeds stability limit.\n";
      }
    if(C2_ode > 1)
      {
        error_message << "C2_ode = " << C2_ode << " exceeds stability limit.\n";
      }
    if(C1_pde > 1)
      {
        error_message << "C1_pde = " << C1_pde << " exceeds stability limit.\n";
      }
    if(C2_pde > 1)
      {
        error_message << "C2_pde = " << C2_pde << " exceeds stability limit.\n";
      }

    if(!error_message.str().empty())
      {
        throw std::runtime_error(error_message.str());
      }
  }

} // namespace stability
} // namespace sim
