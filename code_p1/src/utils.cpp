// By including utils.hpp, we also include all 
// the headers included in utils.hpp
#include "utils.hpp"


int w_file(const std::string& filename, const std::vector<double>& x, const std::vector<double>& u)
{
  // Create an "output file stream" (type std::ofstream)
  // and connect it to our filename.
  std::ofstream ofile;
  ofile.open(filename);

  // Send vectors to output file
  int prec = 3;
  for (size_t i = 0; i < x.size(); i++)
  {
    ofile << std::setprecision(prec) << std::scientific << x[i] << ' ' << u[i] << std::endl;
  }
  
  
  
  // Close the output file
  ofile.close();

  // All is well. Exit program with return code 0.
  return 0;
}


