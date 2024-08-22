// matrix/kronp.cpp
#include "kronp.hpp"
#include <execution> // For parallel execution policies
#include <iostream>
#include <vector>

namespace sim
{
namespace matrix
{

  // Function to calculate the Kronecker product of two matrices
  std::vector<std::vector<double>>
  kron(const std::vector<std::vector<double>> &A,
       const std::vector<std::vector<double>> &B)
  {
    size_t aRows = A.size();
    size_t aCols = A[0].size();
    size_t bRows = B.size();
    size_t bCols = B[0].size();

    std::vector<std::vector<double>> C(aRows * bRows,
                                       std::vector<double>(aCols * bCols, 0));

    // Parallelize the outer loop
    std::for_each(
      std::execution::par_unseq, A.begin(), A.end(), [&](const auto &aRow) {
        size_t i = &aRow - &A[0]; // Get the index of the current row in A
        for(size_t j = 0; j < aCols; ++j)
          {
            for(size_t p = 0; p < bRows; ++p)
              {
                for(size_t q = 0; q < bCols; ++q)
                  {
                    C[i * bRows + p][j * bCols + q] = A[i][j] * B[p][q];
                  }
              }
          }
      });

    return C;
  }

  // Function to calculate the Kronecker power of a matrix
  std::vector<std::vector<double>>
  kronp(const std::vector<std::vector<double>> &A, int k)
  {
    std::vector<std::vector<double>> B = {{1.0}};

    for(int i = 0; i < k; ++i)
      {
        B = kron(B, A);
      }

    return B;
  }

} // namespace matrix
} // namespace sim
