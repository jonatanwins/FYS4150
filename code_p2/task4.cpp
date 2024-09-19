#include "utils.hpp" 


int main() {

    int N = 6;
    double n = N + 1;
    double h = 1/n;
    double a = -1/std::pow(h,2);
    double d = 2/std::pow(h,2);

    arma::mat A = create_tridiagonal(N, a, d);

    double eps = 1e-4;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1e5;
    int iterations = 0;
    bool converged = true;

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
    
    // Getting sorted indices for eigenvalues in decending order (as Armadillo and Analytical set)
    arma::uvec sort_indices = arma::sort_index(eigenvalues, "ascend");

    // Sorting eigenvalues and eigenvectors 
    arma::vec sorted_eigval = eigenvalues(sort_indices);
    arma::mat sorted_eigvec = eigenvectors.cols(sort_indices);

    sorted_eigval.print("Jacobi Rotation Eigvalues");
    sorted_eigvec.print("Jacobi Rotation Eigenvectors");
}

    