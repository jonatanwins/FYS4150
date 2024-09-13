#include "utils.hpp" 

int main()
{
    /*
    int N = 6;
    double n = N + 1;
    double h = 1/n;

    std::cout << "N = " << N << ", n = " << n << ", h = " << h << std::endl;

    double a = -1/std::pow(h,2);
    double d = 2/std::pow(h,2);

    arma::mat A(N, N, arma::fill::zeros); // Matrix (N x N)

    
    // Fill A -------------------------------------------------------------------------
    for (size_t i = 0; i < A.n_rows; i++){
        A(i,i) = d; // Main diagonal
        if (i > 0)
        {
            A(i,i-1) = a; // Lower diagonal
            A(i-1,i) = a; // Upper diagonal
        }
    }
    */

    int N = 99;
    double n = N + 1;
    double h = 1/n;
    double a = -1/std::pow(h,2);
    double d = 2/std::pow(h,2);

    arma::mat A = create_tridiagonal(N, a, d);

    //A.print("Matrix A:");

    // Calculate eigenvals and eigenvec--------------------------------------------------
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A); 
    // Return eigenvectors of unit norm 

    //eigval.print("Eigenvalues:");
    //eigvec.print("Eigenvectors:");

    // Compare with analytical eigenvals and eigenvec -----------------------------------
    arma::vec a_eigval;
    arma::mat a_eigvec;

    std::tie(a_eigval, a_eigvec) = analytical(N, d, a);
    //a_eigval.print("Analytical Eigvalues");
    //a_eigvec.print("Analytical Eigenvectors");

    
    // Write data to output file
    std::string filename = "data/analytical_eigvec_100.txt";
    w_file(filename, a_eigvec);

    return 0;
}

// Run Windows: g++ task2.cpp -o task2.exe -L/ -larmadillo