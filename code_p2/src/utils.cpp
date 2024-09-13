// By including utils.hpp, we also include all 
// the headers included in utils.hpp
#include "utils.hpp"

// Create tridiagonal Matrix
arma::mat create_tridiagonal(int& N,  double& a,  double& d)
{
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

  return A;
}

// Calculating eigenvals and eigenvec analytical
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
    eigvec =  - arma::normalise(eigvec);  

    return {eigval, eigvec};
}

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
    assert(A.is_square());
    assert(A.is_symmetric());

    double max_value = 0;

    for (int i = 0; i < A.n_rows; i++) {
        for (int j=i+1; j < A.n_cols; j++) {
            if (std::abs(A(i,j)) > std::abs(max_value)) {
                max_value = A(i,j);
                k = i;
                l = j;
                // std::cout << "max_value: " << max_value << " k: " << k << " l: " << l << std::endl;
                // std::cout << "A" << A << std::endl;
            }
        }
    }
    return max_value;
}

int w_file(const std::string& filename, const arma::mat& A)
{
    // Create an "output file stream" (type std::ofstream)
    // and connect it to our filename.
    std::ofstream ofile;
    ofile.open(filename);

    //A.print();

    // Send data to output file
    for (size_t i = 0; i < A.n_rows; i++)
    {
        for (size_t j = 0; j < A.n_cols; j++)
        {
            ofile << A(i,j) << " "; 
        }
        
        ofile << std::endl;
        
    }
    
    //ofile << A.print() << std::endl;

    /*
    int prec = 10;
    for (size_t i = 0; i < .size(); i++)
    {
    ofile << std::setprecision(prec) << std::scientific << x[i] << ' ' << u[i] << std::endl;
    } */

    // Close the output file
    ofile.close();

    return 0;
}

