#include "Solver.hpp"
#include <stdexcept>

namespace sim::solvers
{

Eigen::VectorXd
Solver::interp1(const std::vector<double> &ts, const Eigen::MatrixXd &F0,
                double t) const
{
  int n = ts.size();

  if(t < ts.front() || t > ts.back())
    {
      throw std::out_of_range("t is out of bounds");
    }

  int i = 0;
  while(i < n - 1 && t >= ts[i + 1])
    {
      ++i;
    }

  if(t == ts[i])
    {
      return F0.row(i);
    }
  else if(t == ts[i + 1])
    {
      return F0.row(i + 1);
    }

  double t1 = ts[i];
  double t2 = ts[i + 1];
  double ratio = (t - t1) / (t2 - t1);

  return (1 - ratio) * F0.row(i) + ratio * F0.row(i + 1);
}

} // namespace sim::solvers
