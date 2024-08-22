#include "GetPot" // for reading parameters
#include "MainSimulation.hpp"
#ifdef GNUPLOT
#include "gnuplot-iostream.hpp" // interface with gnuplot
#endif
#include "params/SimulationParameters.hpp"
#include "params/readParameters.hpp"
#include <algorithm> // For std::generate
#include <cmath>     // (for sqrt)
#include <iostream>  // input output
#include <numeric>   // For std::iota
#include <tuple>     // I am using tuples to pass data to gnuplot iostream
#include <vector>    // For vectors

/*!
  @file main.cpp
  @brief Carleman Burgers

  @detail
% Solution of the inhomogeneous, viscous Burgers equation, possibly with
% linear damping, using a direct application of the Carleman method
% combined with Euler's method. The result is compared with solutions from
% inbuilt MATLAB solvers.
%
% This script in its current state should reproduce the results in
% https://arxiv.org/abs/2011.03185, "Efficient quantum algorithm for
% dissipative nonlinear differential equations" by Jin-Peng Liu,
% Herman Øie Kolden, Hari K. Krovi, Nuno F. Loureiro, Konstantina Trivisa,
% Andrew M. Childs.
 */

//! a helper function. Prints a help message
void
printHelp()
{
  std::cout
    << "USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"
    << std::endl;
  std::cout << "-h this help" << std::endl;
  std::cout << "-v verbose output" << std::endl;
}

//! main program
int
main(int argc, char **argv)
{
  bool   jsonfile = false;
  int    status{0}; // Status of the program
  GetPot cl(argc, argv);
  if(cl.search(2, "-h", "--help"))
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose = cl.search(1, "-v");
  // Get file with parameter values
  std::string filename = cl.follow("data/parameters.pot", "-p");
  auto        pos = filename.find(".json");
  if(pos != std::string::npos)
    {
      jsonfile = true;
      std::cout << "Json input file\n";
    }
  else
    {
      jsonfile = false;
      std::cout << "Getpot input file\n";
    }
  std::cout << "Reading parameters from " << filename << std::endl;
  sim::params::SimulationParameters param;
  if(jsonfile)
    param = sim::params::readParameters_json(filename, verbose);
  else
    param = sim::params::readParameters(filename, verbose);

  // Pass param as a reference to MainSimulation
  sim::MainSimulation simulation(param);

  simulation.initialize();
  simulation.run();
  return 0;
}
