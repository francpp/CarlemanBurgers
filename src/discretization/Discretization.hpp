#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include "params/SimulationParameters.hpp"
#include <vector>

namespace sim::discretization
{

/**
 * @class Discretization
 * @brief Handles the spatial and temporal discretization of the simulation
 * domain.
 */
class Discretization
{
public:
  /**
   * @brief Constructor for the Discretization class.
   * @param params Reference to the simulation parameters.
   */
  Discretization(const sim::params::SimulationParameters &params);

  /**
   * @brief Creates the spatial and temporal discretization based on the
   * simulation parameters.
   */
  void createDiscretization();

  // Getters for discretized values
  /**
   * @brief Gets the spatial discretization points for the main simulation.
   * @return Reference to a vector containing the spatial discretization points.
   */
  const std::vector<double> &
  getXs() const
  {
    return xs;
  }

  /**
   * @brief Gets the temporal discretization points for the main simulation.
   * @return Reference to a vector containing the temporal discretization
   * points.
   */
  const std::vector<double> &
  getTs() const
  {
    return ts;
  }

  /**
   * @brief Gets the spatial discretization points for the PDE solver.
   * @return Reference to a vector containing the spatial discretization points
   * for the PDE solver.
   */
  const std::vector<double> &
  getXsPde() const
  {
    return xs_pde;
  }

  /**
   * @brief Gets the temporal discretization points for the PDE solver.
   * @return Reference to a vector containing the temporal discretization points
   * for the PDE solver.
   */
  const std::vector<double> &
  getTsPde() const
  {
    return ts_pde;
  }

  /**
   * @brief Gets the temporal discretization points for the ODE solver.
   * @return Reference to a vector containing the temporal discretization points
   * for the ODE solver.
   */
  const std::vector<double> &
  getTsOde() const
  {
    return ts_ode;
  }

  // Getter for SimulationParameters
  /**
   * @brief Gets the reference to the simulation parameters.
   * @return Reference to the SimulationParameters object.
   */
  const sim::params::SimulationParameters &
  getParams() const
  {
    return params;
  }

private:
  const sim::params::SimulationParameters
                     &params; ///< Reference to the simulation parameters.
  std::vector<double> xs, ts, xs_pde, ts_pde,
    ts_ode; ///< Vectors for storing discretization points.

  double dx, dt, dx_pde, dt_pde,
    dt_ode; ///< Step sizes for spatial and temporal discretization.
};

/**
 * @brief Overloads the << operator to print the discretization details.
 * @param out Output stream.
 * @param d Discretization object.
 * @return Reference to the output stream.
 */
std::ostream &operator<<(std::ostream &, const Discretization &);

} // namespace sim::discretization

#endif // DISCRETIZATION_HPP
