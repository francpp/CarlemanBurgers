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
  double x0 = -params.L0 / 2;
  double x1 = params.L0 / 2;
  double t0 = 0;
  double t1 = params.T;

  dx = (x1 - x0) / (params.nx - 1);
  dt = (t1 - t0) / (params.nt - 1);

  xs.resize(params.nx);
  ts.resize(params.nt);

  for(int i = 0; i < params.nx; ++i)
    xs[i] = x0 + i * dx;
  for(int i = 0; i < params.nt; ++i)
    ts[i] = t0 + i * dt;

  int    nt_ode = params.nt * 10;
  double dt_ode = (t1 - t0) / (nt_ode - 1);
  ts_pde.resize(params.nt_pde);
  xs_pde.resize(params.nx_pde);

  for(int i = 0; i < params.nx_pde; ++i)
    xs_pde[i] = x0 + i * dx;
  for(int i = 0; i < params.nt_pde; ++i)
    ts_pde[i] = t0 + i * dt;
}

} // namespace sim::discretization
