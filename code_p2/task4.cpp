#include "utils.hpp" 


int main() {

    int N = 6;
    double n = N + 1;
    double h = 1/n;
    double a = -1/std::pow(h,2);
    double d = 2/std::pow(h,2);

    arma::mat A = create_tridiagonal(N, a, d);

    // A = {
    //     {1, 0, 0, 0.5},
    //     {0, 1, -0.7, 0},
    //     {0, -0.7, 1, 0},
    //     {0.5, 0, 0, 1}
    // };

    double eps = 1e-4;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1e5;
    int iterations = 0;
    bool converged = true;

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
    // std::cout << "Eigenvalues: "<< eigenvalues << std::endl;
    // std::cout << "Eigenvectors: "<< eigenvectors << std::endl;
    // std::cout << "Iterations: "<< iterations << std::endl;
    // std::cout << "Converged: "<< converged << std::endl;

    eigenvalues.print("Jacobi Rotation Eigvalues");
    eigenvectors.print("Jacobi Rotation Eigenvectors");
}
    