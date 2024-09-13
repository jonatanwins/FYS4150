#include "utils.hpp" 

void jacobi_rotate(arma::mat& A, arma::mat& R, int i, int j);

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);