#include "GetPot"
#include "readParameters.hpp"
#include "json.hpp"
#include <fstream>
parameters
readParameters(std::string const &filename, bool verbose)
{
  // Parameter default constructor fills it with the defaults values
  parameters defaults;
  // checks if file exixts and is readable
  std::ifstream check(filename);
  if(!check)
    {
      std::cerr << "ERROR: Parameter file " << filename << " does not exist"
                << std::endl;
      std::cerr << "Reverting to default values." << std::endl;
      if(verbose)
        std::cout << defaults;
      check.close();
      return defaults;
    }
  else
    check.close();

  GetPot     ifile(filename.c_str());
  parameters values;
  // Read parameters from getpot ddata base
  values.nx = ifile("nx", defaults.nx);
  values.nt = ifile("nt", defaults.nt);
  values.nx_pde = ifile("nx_pde", defaults.nx_pde);
  values.nt_pde = ifile("nt_pde", defaults.nt_pde);
  values.Re0 = ifile("Re0", defaults.Re0);
  values.L0 = ifile("L0", defaults.L0);
  values.U0 = ifile("U0", defaults.U0);
  values.beta = ifile("beta", defaults.beta);
  values.f = ifile("f", defaults.f);
  values.T = ifile("T", defaults.T);
  //values.F0_fun = ifile("F0_fun", defaults.F0_fun);
  values.N_max = ifile("N_max", defaults.N_max);
  values.ode_deg = ifile("ode_deg", defaults.ode_deg);
  if(verbose)
    {
      std::cout << "PARAMETER VALUES IN GETPOT FILE"
                << "\n";
      ifile.print();
      std::cout << std::endl;
      std::cout << "ACTUAL VALUES"
                << "\n"
                << values;
    }
  return values;
}

parameters
readParameters_json(std::string const &filename, bool verbose)
{
  // Parameter default constructor fills it with the defaults values
  parameters defaults;
  // checks if file exixts and is readable
  std::ifstream check(filename);
  if(!check)
    {
      std::cerr << "ERROR: Parameter file " << filename << " does not exist"
                << std::endl;
      std::cerr << "Reverting to default values." << std::endl;
      if(verbose)
        std::cout << defaults;
      check.close();
      return defaults;
    }
  else
    check.close();

  std::ifstream jfile(filename);
  nlohmann::json ifile;
  jfile>>ifile;
  parameters values;
  // Read parameters from getpot ddata base
  values.nx = ifile.value("nx", defaults.nx);
  values.nt = ifile.value("nt", defaults.nt);
  values.nx_pde = ifile.value("nx_pde", defaults.nx_pde);
  values.nt_pde = ifile.value("nt_pde", defaults.nt_pde);
  values.Re0 = ifile.value("Re0", defaults.Re0);
  values.L0 = ifile.value("L0", defaults.L0);
  values.U0 = ifile.value("U0", defaults.U0);
  values.beta = ifile.value("beta", defaults.beta);
  values.f = ifile.value("f", defaults.f);
  values.T = ifile.value("T", defaults.T);
  //values.F0_fun = ifile.value("F0_fun", defaults.F0_fun);
  values.N_max = ifile.value("N_max", defaults.N_max);
  values.ode_deg = ifile.value("ode_deg", defaults.ode_deg);
  if(verbose)
    {
      std::cout << "PARAMETER VALUES IN JSON FILE"
                << "\n";
      std::cout<<std::setw(4)<<ifile;
      std::cout << std::endl;
      std::cout << "ACTUAL VALUES"
                << "\n"
                << values;
    }
  return values;
}
