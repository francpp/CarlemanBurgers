#ifndef READ_PARAMETERS_HPP
#define READ_PARAMETERS_HPP

#include "SimulationParameters.hpp"
#include <string>

namespace sim::params
{
//! Reads problem parameters from GetPot or json file
/*!
  @param filename The file with the new values
  @param verbose Prints some information on the parameters
*/
SimulationParameters readParameters(const std::string &filename,
                                    bool               verbose = false);

//! Reads problem parameters from GetPot file
/*!
  @param filename The getopot file with the new values
  @param verbose Prints some information on the parameters
 */
SimulationParameters readParameters_pot(const std::string &filename,
                                        bool               verbose = false);

//! Reads problem parameters from json file
/*!
  @param filename The json file with the parameter values
  @param verbose Prints some information on the parameters
 */
SimulationParameters readParameters_json(const std::string &filename,
                                         bool               verbose = false);

} // namespace sim::params

#endif // READ_PARAMETERS_HPP
