#ifndef INITIAL_CONDITIONS_HPP
#define INITIAL_CONDITIONS_HPP

#include "discretization/Discretization.hpp"
#include "params/SimulationParameters.hpp"
#include <functional>
#include <vector>

namespace sim::initial_conditions
{

class InitialConditions
{
public:
  InitialConditions(const sim::params::SimulationParameters   &params,
                    const sim::discretization::Discretization &discretization);

  void computeInitialConditions();
  const std::vector<double> &
  getU0s() const
  {
    return u0s;
  }

private:
  const sim::params::SimulationParameters   &params;
  const sim::discretization::Discretization &discretization;
  std::vector<double>                        u0s;
};

} // namespace sim::initial_conditions

#endif // INITIAL_CONDITIONS_HPP
