// Now we include headers we need
#include <iostream>
#include <iomanip> // scientific, setprecision()
#include <vector>
#include <string>
#include <armadillo> 
#include <cmath> // For M_PI



// Function *declarations*.


int main()
{
    int N = 6;
    double n = N + 1;
    double h = 1/n;

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

    A.print("Matrix A:");
    
    std::cout << "N = " << N << ", n = " << n << ", h = " << h << std::endl;

    // Calculate eigenvals and eigenvec--------------------------------------------------
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A); 
    // Return eigenvectors of unit norm 

    eigval.print("Eigenvalues:");
    eigvec.print("Eigenvectors:");

    // Compare with analytical eigenvals and eigenvec -----------------------------------
    arma::vec a_eigval;
    arma::mat a_eigvec;

    std::tie(a_eigval, a_eigvec) = analytical(N, d, a);
    a_eigval.print("Analytical Eigvalues");
    a_eigvec.print("Analytical Eigenvectors");
    return 0;
}


std::pair<arma::vec, arma::mat> analytical(int N, double d, double a)
{
    arma::vec eigval(N);
    arma::mat eigvec(N, N);

    for (size_t j = 0; j < N; j++)
    {
        eigval(j) = d + 2*a*cos((j+1) * M_PI / (N + 1));
        for (size_t i = 0; i < N; i++)
        {
            eigvec(i,j) = sin((i+1) * (j+1) * M_PI / (N + 1));          
        }
        
    } 
    // Normalize unit norm, use negative values to fit with armadillo values
    eigvec = - arma::normalise(eigvec);  

    return {eigval, eigvec};
}

// Run: g++ task2.cpp -o task2.exe -L/ -larmadillo