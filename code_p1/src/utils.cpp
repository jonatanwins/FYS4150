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
  int prec = 10;
  for (size_t i = 0; i < x.size(); i++)
  {
    ofile << std::setprecision(prec) << std::scientific << x[i] << ' ' << u[i] << std::endl;
  } 
  
  // Close the output file
  ofile.close();

  return 0;
}

int w_file_one(const std::string& filename, const std::vector<int>& x, const std::vector<double>& u, const std::vector<double>& v)
{
  std::ofstream ofile;
  ofile.open(filename);

  int prec = 10;
  for (size_t i = 0; i < x.size(); i++)
  {
    ofile << std::setprecision(prec) << std::scientific << x[i] << ' ' << u[i]  << ' ' << v[i] << std::endl;
  } 
  
  ofile.close();

  return 0;
}

// Forward substitution function

std::pair<std::vector<double>, std::vector<double>> forwardSubstitution(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1 ) {

    int n = b.size()-1;
    double h = 1./n;

    std::vector<double> g(n+2);

    g[1] = h*h*f[1] + v_0;
    g[n] = h*h*f[n] + v_n_plus_1;

    for (int i = 2; i <= n-1; i++) {
        g[i] = h*h*f[i];
    }

    std::vector<double> btilde(n+2);
    std::vector<double> gtilde(n+2);

    btilde[1] = b[1];
    gtilde[1] = g[1];
 
    for (int i = 2; i <= n; i++) {
        double w = (a[i]/btilde[i-1]);
        btilde[i] = b[i] - w*c[i-1];
        gtilde[i] = g[i] - w*gtilde[i-1];
    }

    return {gtilde, btilde};   
}


// Forward substitution function for the special algorithm

std::pair<std::vector<double>, std::vector<double>> forwardSubstitution_special(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1 ) {

    int n = b.size()-1;
    double h = 1./n;

    std::vector<double> g(n+2);

    g[1] = h*h*f[1] + v_0;
    g[n] = h*h*f[n] + v_n_plus_1;

    for (int i = 2; i <= n-1; i++) {
        g[i] = h*h*f[i];
    }

    std::vector<double> btilde(n+2);
    std::vector<double> gtilde(n+2);

    btilde[1] = b[1];
    gtilde[1] = g[1];
 
    for (int i = 2; i <= n; i++) {
        double w = (a[i]/btilde[i-1]);
        btilde[i] = i/(i+1);
        gtilde[i] = g[i] - w*gtilde[i-1];
    }

    return {gtilde, btilde};   
}

// Backward substitution function

std::vector<double> backwardSubstitution(const std::vector<double>& gtilde, const std::vector<double>& btilde, const std::vector<double>& c, double v_0, double v_n_plus_1) {

    int n = c.size();
    std::vector<double> v(n+2);

    v[0] = v_0;
    v[n+1] = v_n_plus_1;

    v[n] = gtilde[n]/btilde[n];

    for (int i = n-1; i >= 1; i--) {
        v[i] = (gtilde[i] - v[i+1]*c[i])/btilde[i];
    }

    return v;
}

// Function to solve the tridiagonal matrix equation

std::vector<double> solveTridiagonalMatrixEquation(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1) {

    std::pair<std::vector<double>, std::vector<double>> forw = forwardSubstitution(a, b, c, f, v_0, v_n_plus_1);
    std::vector<double> gtilde = forw.first;
    std::vector<double> btilde = forw.second;
    return backwardSubstitution(gtilde, btilde, c, v_0, v_n_plus_1);
}

// Exact solution : u(x) = 1 - (1-e^{-10})x - e^{-10x}

std::vector<double> exactSolution(double n) {
    
    std::vector<double> u(n+2);
    for (int i = 0; i <= n+1; i++) {
        u[i] = 1 - (1 - exp(-10))*(i/n) - exp(-10*(i/n));
    }
    return u;
}

// Function to solve the tridiagonal matrix equation for the special case

std::vector<double> solveTridiagonalMatrixEquation_special(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1) {

    std::pair<std::vector<double>, std::vector<double>> forw = forwardSubstitution_special(a, b, c, f, v_0, v_n_plus_1);
    std::vector<double> gtilde = forw.first;
    std::vector<double> btilde = forw.second;
    return backwardSubstitution(gtilde, btilde, c, v_0, v_n_plus_1);
}

