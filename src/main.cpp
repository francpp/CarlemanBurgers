#include "GetPot"               // for reading parameters
#ifdef GNUPLOT
#include "gnuplot-iostream.hpp" // interface with gnuplot
#endif
#include "readParameters.hpp"
#include <cmath>    // (for sqrt)
#include <iostream> // input output
#include <tuple> // I am using tuples to pass data to gnuplot iostream
#include <vector> // For vectors
#include <numeric> // For std::iota
#include <algorithm> // For std::generate
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
  using namespace std; // avoid std::
  int    status{0};    // final program status
  bool jsonfile=false;
  GetPot cl(argc, argv);
  if(cl.search(2, "-h", "--help"))
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose = cl.search(1, "-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot", "-p");
  auto pos = filename.find(".json");
  if(pos != std::string::npos)
  {
    jsonfile=true;
    std::cout<<"Json input file\n";
  }
  else
  {
    jsonfile=false;
    std::cout<<"Getpot input file\n";
  }
  cout << "Reading parameters from " << filename << std::endl;
  parameters param;
    if(jsonfile)
      param=   readParameters_json(filename, verbose);
    else
      param=   readParameters(filename, verbose);

  #if __cplusplus < 201703L
    // This version is perfectly fine and
    // works also if you compile with C++17, but with C++17 you
    // may make things simpler
    // Transfer parameters to local variables, to avoid having towrite every time
    // param.xx. I use references to save memory (not really an issue here, it is
    // just to show a possible  use of references)
    const int nx = 16; // Spatial discretization for Euler's method
    const int nt = 4000; // Temporal discretization for Euler's method
    const int nx_pde = 100; // Spatial discretization for the pdepe solver
    const int nt_pde = 40000; // Temporal discretization for the pdepe solver

    const double Re0 = 20; // Desired Reynolds number
    const double L0 = 1; // Domain length
    const double U0 = 1/sqrt(nx-1); // Initial maximum velocity
    const double beta = 0; // Linear damping coefficient
    const double f = 1; // Number of oscillations of the sinusoidal initial condition inside the domain
    const double T = 3; // Simulation time
    //auto F0_fun = [](double t, double x) { return U0*exp(-(x-L0/4)*(x-L0/4)/(2*(L0/32)*(L0/32)))*cos(2*M_PI*t); }; // Source function.
    const int N_max = 4; // Maximum Carleman truncation level
    const int ode_deg = 2; // Degree of the Carleman ODE, should not be changed

  #else
    // C++17 onwards version. This version works only with at least C++17
    // A oneliner! This is called structured bindings. It works because parameter
    // class is an aggregate!

      const auto &[nx, nt, nx_pde, nt_pde, Re0, L0, U0, beta, f, T, N_max, ode_deg] = param;
      
  #endif

  // Truncation levels
  std::vector<int> Ns(N_max);
  std::iota(Ns.begin(), Ns.end(), 1); 

  const auto nu = U0 * L0 / Re0; // Viscosity
  const auto Tnl = L0 / U0; // Nonlinear time
  const auto t_plot = Tnl / 3; // Time to plot solution

  // Spatial domain edges
  const auto x0 = -L0 / 2;
  const auto x1 = L0 / 2;

  // Temporal domain edges
  const auto t0 = 0;
  const auto t1 = T;

  // Euler's method discretization interval sizes and domains
  const auto dx = (x1 - x0) / (nx - 1);
  const auto dt = (t1 - t0) / (nt - 1);

  std::vector<double> xs(nx);
  std::generate(xs.begin(), xs.end(), [n = 0, x0, dx]() mutable {
      return x0 + n++ * dx;
  });

  std::vector<double> ts(nt);
  std::generate(ts.begin(), ts.end(), [n = 0, t0, dt]() mutable {
      return t0 + n++ * dt;
  });

  // ode45 discretization
  const auto nt_ode = nt * 10; // Make it more accurate than the Euler solution
  const auto dt_ode = (t1 - t0) / (nt_ode - 1);
  std::vector<double> ts_ode(nt_ode);
  std::generate(ts_ode.begin(), ts_ode.end(), [n = 0, t0, dt_ode]() mutable {
      return t0 + n++ * dt_ode;
  });

  // pdepe discretization interval sizes
  const auto dx_pde = (x1 - x0) / (nx_pde - 1); // Spatial discretization interval size for pdepe solver
  const auto dt_pde = (t1 - t0) / (nt_pde - 1);
  const auto xs_pde = [&] {
      std::vector<double> v(nx_pde);
      std::generate(v.begin(), v.end(), [n = 0, x0, dx_pde]() mutable {
          return x0 + n++ * dx_pde;
      });
      return v;
  }();

  const auto ts_pde = [&] {
      std::vector<double> v(nt_pde);
      std::generate(v.begin(), v.end(), [n = 0, t0, dt_pde]() mutable {
          return t0 + n++ * dt_pde;
      });
      return v;
  }();

  // Analytic solution
  std::vector<double> theta(10, 1.0);

  std::vector<double> thetaa(10, 2.0);

  std::vector<double> coor(10);

  cout << "Result file: result.dat" << endl;
  ofstream f_stream("result.dat");
  // In gnuplot lines beginning with # are comments
  // \t writes a tab
  f_stream << "#node coordinate\tcomputed solution\texact solution" << std::endl;
  for(int m = 0; m <= 10; m++)
    {
      f_stream.setf(std::ios::left, std::ios::adjustfield);
      f_stream.width(18);
      f_stream.precision(15);
      f_stream << m << "\t\t" << theta[m] << "\t\t" << thetaa[m] << "\n";
      coor[m] = m;
    }
  // If you have gnuplot iostream you get the plot on the screen
#ifdef GNUPLOT
  Gnuplot gp; // gnuplot iostream! Plots solution on the screen
  // It may not work on virtual machines. Take it out in that case
  // Using temporary files (another nice use of tie)
  // Comment this statement if you are not using gnuplot iostream
  // to plot the solution directly on the terminal
  gp << "plot" << gp.file1d(std::tie(coor, theta)) << "w lp lw 2 title 'uh',"
     << gp.file1d(std::tie(coor, thetaa)) << "w l lw 2 title 'uex'"
     << std::endl;
#endif
  f_stream.close();
  return status;
}
