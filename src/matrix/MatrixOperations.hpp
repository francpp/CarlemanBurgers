#ifndef MATRIX_OPERATIONS_HPP
#define MATRIX_OPERATIONS_HPP

#include "discretization/Discretization.hpp"
#include "params/SimulationParameters.hpp"
#include <vector>

namespace sim::matrix
{
class MatrixOperations
{
public:
  std::vector<std::vector<double>> F0, F1, F2;

  MatrixOperations(const params::SimulationParameters   &params,
                   const discretization::Discretization &discretization);
  void createMatrices();
};
} // namespace sim::matrix

#endif // MATRIX_OPERATIONS_HPP
