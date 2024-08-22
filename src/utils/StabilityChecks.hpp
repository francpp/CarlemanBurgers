#ifndef STABILITY_CHECKS_HPP
#define STABILITY_CHECKS_HPP

namespace sim
{
namespace stability
{

  // Function to check CFL conditions and throw an error if any condition is
  // violated.
  void checkCFLConditions(double U0, double dt, double dx, double nu,
                          double dt_ode, double dx_pde, double dt_pde);

} // namespace stability
} // namespace sim

#endif // STABILITY_CHECKS_HPP
