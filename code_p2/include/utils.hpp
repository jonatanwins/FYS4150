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
#include <armadillo>
#include <assert.h>



// Function *declarations*.
std::pair<arma::vec, arma::mat> analytical(int N, double d, double a);
arma::mat create_tridiagonal(int& N,  double& a,  double& d);
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);
void jacobi_rotate(arma::mat& A, arma::mat& R, int i, int j);
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged);
int w_file(const std::string& filename, const arma::mat& A);

#endif  // end of include guard __utils_hpp__