#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include "params/SimulationParameters.hpp"
#include <vector>

namespace sim::discretization
{

class Discretization
{
public:
  Discretization(const sim::params::SimulationParameters &params);

  void createDiscretization();

  // Getters for discretized values
  const std::vector<double> &
  getXs() const
  {
    return xs;
  }
  const std::vector<double> &
  getTs() const
  {
    return ts;
  }
  const std::vector<double> &
  getXsPde() const
  {
    return xs_pde;
  }
  const std::vector<double> &
  getTsPde() const
  {
    return ts_pde;
  }

private:
  const sim::params::SimulationParameters &params;
  std::vector<double>                      xs, ts, xs_pde, ts_pde;

  double dx, dt, dx_pde, dt_pde;
};

} // namespace sim::discretization

#endif // DISCRETIZATION_HPP
