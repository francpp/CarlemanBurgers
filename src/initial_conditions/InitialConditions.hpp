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

  void computeInitialConditions();         // Computes the initial conditions
  void computeForcingBoundaryConditions(); // Computes F0, F1, F2 matrices

  const std::vector<double> &
  getU0s() const
  {
    return u0s;
  }
  const std::vector<std::vector<double>> &
  getF0() const
  {
    return F0;
  }
  const std::vector<std::vector<double>> &
  getF1() const
  {
    return F1;
  }
  const std::vector<std::vector<double>> &
  getF2() const
  {
    return F2;
  }

private:
  const sim::params::SimulationParameters   &params;
  const sim::discretization::Discretization &discretization;

  std::vector<double>              u0s;        // Initial condition array
  std::vector<std::vector<double>> F0, F1, F2; // Forcing matrices

  std::function<double(double, double)> F0_fun; // Source function
};

} // namespace sim::initial_conditions

#endif // INITIAL_CONDITIONS_HPP
