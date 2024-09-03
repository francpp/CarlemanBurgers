#ifndef STABILITY_CHECKS_HPP
#define STABILITY_CHECKS_HPP

namespace sim
{
namespace stability
{

  /**
   * @brief Checks the CFL (Courant–Friedrichs–Lewy) conditions for the
   * simulation.
   *
   * This function checks whether the CFL conditions are satisfied for both
   * explicit Euler, ODE, and PDE discretizations. If any condition is violated,
   * it throws a runtime error.
   *
   * @param U0 Initial velocity.
   * @param dt Time step for the explicit Euler method.
   * @param dx Spatial step size for the explicit Euler method.
   * @param nu Kinematic viscosity.
   * @param dt_ode Time step for the ODE solver.
   * @param dx_pde Spatial step size for the PDE solver.
   * @param dt_pde Time step for the PDE solver.
   * @throws std::runtime_error if any CFL condition is violated.
   */
  void checkCFLConditions(double U0, double dt, double dx, double nu,
                          double dt_ode, double dx_pde, double dt_pde);

} // namespace stability
} // namespace sim

#endif // STABILITY_CHECKS_HPP
