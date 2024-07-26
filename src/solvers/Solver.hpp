#ifndef SOLVER_HPP
#define SOLVER_HPP

namespace sim::solvers
{
class Solver
{
public:
  virtual void solve() = 0;
};
} // namespace sim::solvers

#endif // SOLVER_HPP
