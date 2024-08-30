#include <fstream>
#include <iostream>
#include <vector>

void
plotGraph(const std::vector<double> &x, const std::vector<double> &y)
{
  std::ofstream temp("temp.dat");
  for(size_t i = 0; i < x.size(); ++i)
    {
      temp << x[i] << " " << y[i] << std::endl;
    }
  temp.close();

  // Command to run gnuplot and save the output to a PNG file
  std::string command = "gnuplot -e \"set terminal png; set output "
                        "'output.png'; plot 'temp.dat' with lines\"";
  system(command.c_str());
}

int
main()
{
  std::vector<double> x = {1, 2, 3, 4, 5};
  std::vector<double> y = {1, 4, 9, 16, 25};

  plotGraph(x, y);

  return 0;
}
