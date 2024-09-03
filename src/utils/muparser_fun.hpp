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
  /**
   * @class MuparserFun
   * @brief Wrapper around the muParser library to simplify function parsing and
   * evaluation.
   *
   * The MuparserFun class allows you to define mathematical expressions with
   * variables, and then evaluate those expressions by providing values for the
   * variables.
   */
  class MuparserFun
  {
  public:
    /**
     * @brief Copy constructor for MuparserFun.
     * @param m The MuparserFun object to copy from.
     */
    MuparserFun(const MuparserFun &m) : m_parser(m.m_parser), m_vars(m.m_vars)
    {
      for(auto &var : m_vars)
        {
          m_parser.DefineVar(var.first, var.second);
        }
    }

    /**
     * @brief Constructor that initializes the muParser object with a given
     * expression and variables.
     * @param expr The mathematical expression to parse.
     * @param variables A map of variable names to pointers to their values.
     */
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

    /**
     * @brief Evaluates the expression with the provided values for variables.
     * @param values A map of variable names to their corresponding values to be
     * used in the expression.
     * @return The result of evaluating the expression.
     * @throws mu::Parser::exception_type if the evaluation fails.
     */
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
    std::map<std::string, double *>
      m_vars; ///< Map of variable names to their corresponding pointers.
    mu::Parser
      m_parser; ///< muParser object for parsing and evaluating expressions.
  };
} // namespace utils
} // namespace sim

#endif // MUPARSER_FUN_HPP
