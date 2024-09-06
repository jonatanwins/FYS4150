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



// Below we give some function *declarations*.
// The function *definitions* (the actual code) 
// lives in src/utils.cpp


int w_file(const std::string& filename, const std::vector<double>& x, const std::vector<double>& u);

#endif  // end of include guard __utils_hpp__