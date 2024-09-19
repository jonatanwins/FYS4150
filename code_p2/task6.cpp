#include "utils.hpp" 

int main()
{
    // for n = 10 set N = 9 and change filename and filename2
    // for n = 100 set N = 99 and change filename and filename2
    int N = 9; 
    double n = N + 1;
    double h = 1/n;
    double a = -1/std::pow(h,2);
    double d = 2/std::pow(h,2);

    arma::mat A = create_tridiagonal(N, a, d);

    // Find Analytical eigenvals and eigenvec  -----------------------------------
    arma::vec a_eigval;
    arma::mat a_eigvec;

    std::tie(a_eigval, a_eigvec) = analytical(N, d, a);
    
    // Write data to output file to 
    std::string filename = "data/analytical_eigvec_10.txt";
    w_file(filename, a_eigvec);




    // Find Jacobi eigenvals and eigenvec     -----------------------------------
    arma::vec j_eigval;
    arma::mat j_eigvec;

    double eps = 1e-3; 
    int maxiter = 1e5;
    int iterations = 0;
    bool converged = true;

    jacobi_eigensolver(A, eps, j_eigval, j_eigvec, maxiter, iterations, converged);

    // Getting sorted indices for eigenvalues in decending order
    arma::uvec sort_indices = arma::sort_index(j_eigval, "ascend");

    // Sorting eigenvalues and eigenvectors 
    arma::vec sorted_eigval = j_eigval(sort_indices);
    arma::mat sorted_eigvec = j_eigvec.cols(sort_indices);

    
    // Write data to output file to 
    std::string filename2 = "data/jacobi_eigvec_10.txt";
    w_file(filename2, a_eigvec);


    return 0;
}

// Run Windows: g++ task2.cpp -o task2.exe -L/ -larmadillo