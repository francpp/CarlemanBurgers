#ifndef MUPARSER_FUN_HPP
#define MUPARSER_FUN_HPP

#include <muParser.h>
#include <iostream>
#include <map>
#include <string>

namespace sim
{
namespace utils
{
  class MuparserFun
  {
  public:
    MuparserFun(const MuparserFun &m) : m_parser(m.m_parser), m_vars(m.m_vars)
    {
      for(auto &var : m_vars)
        {
          m_parser.DefineVar(var.first, var.second);
        }
    }

    MuparserFun(const std::string                     &expr,
                const std::map<std::string, double *> &variables)
    {
      try
        {
          for(const auto &var : variables)
            {
              m_parser.DefineVar(var.first, var.second);
              m_vars[var.first] = var.second;
            }
          m_parser.DefineConst("pi", 3.141592653589793); // Define pi
          m_parser.SetExpr(expr);
        }
      catch(mu::Parser::exception_type &e)
        {
          std::cerr << "Error in MuparserFun constructor: " << e.GetMsg()
                    << std::endl;
          throw; // Re-throw to handle higher up if needed
        }
    }

    double
    operator()(const std::map<std::string, double> &values)
    {
      try
        {
          for(const auto &val : values)
            {
              *m_vars[val.first] = val.second;
            }
          return m_parser.Eval();
        }
      catch(mu::Parser::exception_type &e)
        {
          std::cerr << "Error in MuparserFun operator(): " << e.GetMsg()
                    << std::endl;
          throw; // Re-throw to handle higher up if needed
        }
    }

  private:
    std::map<std::string, double *> m_vars;
    mu::Parser                      m_parser;
  };
} // namespace utils
} // namespace sim

#endif // MUPARSER_FUN_HPP
