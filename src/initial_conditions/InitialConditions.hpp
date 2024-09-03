#ifndef INITIAL_CONDITIONS_HPP
#define INITIAL_CONDITIONS_HPP

#include "discretization/Discretization.hpp"
#include "params/SimulationParameters.hpp"
#include <functional>
#include <vector>

namespace sim::initial_conditions
{

/**
 * @class InitialConditions
 * @brief Class responsible for computing the initial conditions and forcing
 * boundary conditions for the simulation.
 */
class InitialConditions
{
public:
  /**
   * @brief Constructor for InitialConditions.
   * @param params Reference to the simulation parameters.
   * @param discretization Reference to the discretization object.
   */
  InitialConditions(const sim::params::SimulationParameters   &params,
                    const sim::discretization::Discretization &discretization);

  /**
   * @brief Computes the initial conditions based on the given expressions.
   */
  void computeInitialConditions();

  /**
   * @brief Computes the forcing boundary conditions for the simulation.
   *        This includes the calculation of matrices F0, F1, and F2.
   */
  void computeForcingBoundaryConditions();

  /**
   * @brief Computes the matrix F1 which contains the linear coefficients
   * corresponding to the linear parts of the equation, such as the diffusion
   * term.
   */
  void computeF1();

  /**
   * @brief Computes the matrix F2 which captures the quadratic interactions
   * between the components of u, representing contributions from the advection
   * term.
   */
  void computeF2();

  /**
   * @brief Gets the initial condition values.
   * @return Reference to the vector containing initial condition values.
   */
  const std::vector<double> &
  getU0s() const
  {
    return u0s;
  }

  /**
   * @brief Gets the forcing matrix F0.
   * @return Reference to the 2D vector containing F0 matrix values.
   */
  const std::vector<std::vector<double>> &
  getF0() const
  {
    return F0;
  }

  /**
   * @brief Gets the matrix F1.
   * @return Reference to the 2D vector containing F1 matrix values.
   */
  const std::vector<std::vector<double>> &
  getF1() const
  {
    return F1;
  }

  /**
   * @brief Gets the matrix F2.
   * @return Reference to the 2D vector containing F2 matrix values.
   */
  const std::vector<std::vector<double>> &
  getF2() const
  {
    return F2;
  }

private:
  const sim::params::SimulationParameters
    &params; ///< Reference to the simulation parameters.
  const sim::discretization::Discretization
    &discretization; ///< Reference to the discretization object.

  std::vector<double>              u0s;        ///< Initial condition array.
  std::vector<std::vector<double>> F0, F1, F2; ///< Forcing matrices.

  std::string F0_expr; ///< Source function expression.
  std::string U0_expr; ///< Initial condition function expression.
};

/**
 * @brief Overload of the << operator to print the initial conditions and
 * matrices.
 * @param out Output stream.
 * @param d InitialConditions object.
 * @return Reference to the output stream.
 */
std::ostream &operator<<(std::ostream &, const InitialConditions &);

} // namespace sim::initial_conditions

#endif // INITIAL_CONDITIONS_HPP
