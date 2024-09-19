#include "utils.hpp" 


int main() {
    std::vector<int> N = {5, 10, 15, 20, 30, 40, 50, 100, 150};
    double eps = 1e-3;
    arma::vec eigenvalues;
    arma::mat eigenvectors;

    int maxiter = 1e5;

    for (int i = 0; i < N.size(); i++) {
        double n = N[i] + 1;
        double h = 1/n;
        double a = -1/std::pow(h,2);
        double d = 2/std::pow(h,2);

        // arma::mat A = create_tridiagonal(N[i], a, d);

        // Generate random N*N matrix
        arma::mat A = arma::mat(N[i], N[i]).randn();  

        // Symmetrize the matrix by reflecting the upper triangle to lower triangle
        A = arma::symmatu(A);  

        int iterations = 0;
        bool converged = true;

        // A.print();
        jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

        std::cout << "N: "<< N[i] << " Iterations: "<< iterations << std::endl;
        std::cout << "Converged: "<< converged << std::endl;
    }
}