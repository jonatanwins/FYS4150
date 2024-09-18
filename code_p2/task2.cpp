#include "utils.hpp" 

int main()
{
    int N = 6;
    double n = N + 1;
    double h = 1/n;
    double a = -1/std::pow(h,2);
    double d = 2/std::pow(h,2);

    arma::mat A = create_tridiagonal(N, a, d);
    A = {
        {1, 0, 0, 0.5},
        {0, 1, -0.7, 0},
        {0, -0.7, 1, 0},
        {0.5, 0, 0, 1}
    };

    // Calculate eigenvals and eigenvec form Armadillo --------------------------------------------------
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A); 
    // Return eigenvalues eigenvectors of unit norm from armadillo
    eigval.print("Armadillo Eigvalues");
    eigvec.print("Armadillo Eigenvectors");


    // Compare with analytical eigenvals and eigenvec -----------------------------------
    arma::vec a_eigval;
    arma::mat a_eigvec;

    std::tie(a_eigval, a_eigvec) = analytical(N, d, a);
    // Printing analytical eigenvalues and eigenvectors
    a_eigval.print("Analytical Eigvalues");
    a_eigvec.print("Analytical Eigenvectors");

    return 0;
}

// Run Windows: g++ task2.cpp -o task2.exe -L/ -larmadillo