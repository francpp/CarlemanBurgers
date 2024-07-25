#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <cmath> // for sqrt
#include <iosfwd>
/*!
 * A structure holding the parameters
 *
 * It is an aggregate, you can use structured binding and brace initialization
 */
struct parameters
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
  //! Linear damping coefficient
  double beta = 0;
  //! Number of oscillations of the sinusoidal initial condition inside the
  //! domain
  double f = 1;
  //! Simulation time
  double T = 3;
  //! Source function
  // TODO F0_fun
  //! Maximum Carleman truncation level
  int N_max = 4;
  //! Degree of the Carleman ODE
  int ode_deg = 2;
};
//! Prints parameters
std::ostream &operator<<(std::ostream &, const parameters &);
#endif
