// First we add an "include guard". It ensures that this header file can 
// only end up being included *once* for the compilation of any given .cpp file
#ifndef __utils_hpp__  
#define __utils_hpp__

// Now we include headers we need
#include <iostream>
#include <iomanip> // scientific, setprecision()
#include <vector>
#include <string>
#include <cmath>   // exp()
#include <fstream> // ofstream   (read text)
#include <chrono>



// Below we give some function *declarations*.
// The function *definitions* (the actual code) 
// lives in src/utils.cpp


int w_file(const std::string& filename, const std::vector<double>& x, const std::vector<double>& u);
int w_file_one(const std::string& filename, const std::vector<int>& x, const std::vector<double>& u);

std::pair<std::vector<double>, std::vector<double>> forwardSubstitution(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1 );
std::vector<double> backwardSubstitution(const std::vector<double>& gtilde, const std::vector<double>& btilde, const std::vector<double>& c, double v_0, double v_n_plus_1);
std::vector<double> solveTridiagonalMatrixEquation(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1);
std::vector<double> exactSolution(double n);


#endif  // end of include guard __utils_hpp__