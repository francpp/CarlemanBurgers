// matrix/kronp.hpp
#ifndef KRONP_HPP
#define KRONP_HPP

#include <vector>

namespace sim
{
namespace matrix
{

  // Function to calculate the Kronecker product of two matrices
  std::vector<std::vector<double>>
  kron(const std::vector<std::vector<double>> &A,
       const std::vector<std::vector<double>> &B);

  // Function to calculate the Kronecker power of a matrix
  std::vector<std::vector<double>>
  kronp(const std::vector<std::vector<double>> &A, int k);

} // namespace matrix
} // namespace sim

#endif // KRONP_HPP
