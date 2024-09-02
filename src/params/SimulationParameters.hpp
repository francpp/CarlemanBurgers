#ifndef SIMULATION_PARAMETERS_HPP
#define SIMULATION_PARAMETERS_HPP

#include <cmath> // for sqrt
#include <functional>
#include <iosfwd>
#include <iostream>
#include <string>

namespace sim::params
{

/*!
 * A structure holding the parameters
 *
 * It is an aggregate, you can use structured binding and brace initialization
 */
struct SimulationParameters
{
  //! Spatial discretization for Euler's method
  int nx = 16;
  //! Temporal discretization for Euler's method
  int nt = 4000;
  //! Spatial discretization for the pdepe solver
  int nx_pde = 100;
  //! Temporal discretization for the pdepe solver
  int nt_pde = 40000;
  //! Desired Reynolds number
  double Re0 = 20;
  //! Domain length
  double L0 = 1;
  //! Initial maximum velocity
  double U0 = 1 / sqrt(nx - 1);
  //! Viscosity
  double nu = U0 * L0 / Re0;
  //! Linear damping coefficient
  double beta = 0;
  //! Number of oscillations of the sinusoidal initial condition inside the
  //! domain
  double f = 1;
  //! Simulation time
  double T = 3;
  //! Maximum Carleman truncation level
  int N_max = 4;
  //! Degree of the Carleman ODE
  int ode_deg = 2;
  //! Initial condition function
  std::string U0_fun = "-U0*sin(2*pi*f*x/L0)";
  //! Source function
  std::string F0_fun = "U0*exp(-(x-L0/4)^2/(2*(L0/32)^2))*cos(2*pi*t)";

  void initialize();
};

//! Prints parameters
std::ostream &operator<<(std::ostream &, const SimulationParameters &);

} // namespace sim::params

#endif // SIMULATION_PARAMETERS_HPP
