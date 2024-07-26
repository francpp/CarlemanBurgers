#include "InitialConditions.hpp"
#include <algorithm>
#include <cmath>

namespace sim::initial_conditions
{

InitialConditions::InitialConditions(
  const sim::params::SimulationParameters   &params,
  const sim::discretization::Discretization &discretization)
  : params(params), discretization(discretization)
{}

void
InitialConditions::computeInitialConditions()
{
  const auto &xs = discretization.getXs();
  u0s.resize(xs.size());

  std::transform(xs.begin(), xs.end(), u0s.begin(), [&](double x) {
    return -params.U0 * std::sin(2 * M_PI * params.f * x / params.L0);
  });
}

} // namespace sim::initial_conditions
