#include "utils.hpp" 


int main() {
    std::vector<int> N = {2, 3, 4, 5, 6, 7, 8, 9, 10};
    double eps = 1e-4;
    arma::vec eigenvalues;
    arma::mat eigenvectors;

    int maxiter = 1e5;

    for (int i = 0; i < N.size(); i++) {
        double n = N[i] + 1;
        double h = 1/n;
        double a = -1/std::pow(h,2);
        double d = 2/std::pow(h,2);

        arma::mat A = create_tridiagonal(N[i], a, d);

        // Generate random N*N matrix
        // arma::mat A = arma::mat(N, N).randn();  

        // Symmetrize the matrix by reflecting the upper triangle to lower triangle
        // A = arma::symmatu(A);  

        int iterations = 0;
        bool converged = true;

        A.print();
        jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

        std::cout << "N: "<< N[i] << std::endl;
        std::cout << "Iterations: "<< iterations << std::endl;
        std::cout << "Converged: "<< converged << std::endl;
    }
}