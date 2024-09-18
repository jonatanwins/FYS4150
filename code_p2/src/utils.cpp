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


void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                const int maxiter, int& iterations, bool& converged)
{
    arma::mat A_m = A;
    arma::mat A_m_plus_1 = A_m;

    int N = A_m.n_rows;

    arma::mat R_m = arma::mat(N, N, arma::fill::eye);
    arma::mat R_m_plus_1 = R_m;

    int k = 0;
    int l = 0;
    max_offdiag_symmetric(A_m, k, l); // k and l are updated by reference

    while (std::abs(A_m(k,l)) > eps)
    {
        if (iterations > maxiter)
        {
            converged = false;
            break;
        }

        double tau = (A_m(l,l) - A_m(k,k)) / (2*A_m(k,l)) ;
        double t_1 = 1/(tau + std::sqrt(1 + std::pow(tau,2)));
        double t_2 = 1/(tau - std::sqrt(1 + std::pow(tau,2)));
        double t;
        if (tau < 0)
        {
            t = t_1;
        }
        else
        {
            t = t_2;
        }

        double c = 1/std::sqrt(1 + std::pow(t,2));
        double s = c*t;

        A_m_plus_1(k,k) = A_m(k,k)*std::pow(c,2) - 2*A_m(k,l)*c*s + A_m(l,l)*std::pow(s,2);
        A_m_plus_1(l,l) = A_m(l,l)*std::pow(c,2) + 2*A_m(k,l)*c*s + A_m(k,k)*std::pow(s,2);

        for (int i = 0; i < N; i++) {
            if (i != k && i != l) {
                A_m_plus_1(i,k) = A_m(i,k)*c - A_m(i,l)*s;
                A_m_plus_1(k,i) = A_m_plus_1(i,k);
                A_m_plus_1(i,l) = A_m(i,l)*c + A_m(i,k)*s;
                A_m_plus_1(l,i) = A_m_plus_1(i,l);
            } 
        }

        A_m_plus_1(k,l) = 0;
        A_m_plus_1(l,k) = 0;
        
        A_m = A_m_plus_1;

        for (int i = 0; i < N; i++) {
                R_m_plus_1(i,k) = R_m(i,k)*c - R_m(i,l)*s;
                R_m_plus_1(i,l) = R_m(i,l)*c + R_m(i,k)*s;
            }
        
        R_m = R_m_plus_1;
        iterations = iterations + 1;

        double max_value = max_offdiag_symmetric(A_m, k, l);

    }

    eigenvalues = arma::diagvec(A_m);
    eigenvectors = R_m; 

}



int w_file(const std::string& filename, const arma::mat& A)
{
    // Create an "output file stream" (type std::ofstream)
    // and connect it to our filename.
    std::ofstream ofile;
    ofile.open(filename);

    // Send data to output file
    for (size_t i = 0; i < A.n_rows; i++)
    {
        for (size_t j = 0; j < A.n_cols; j++)
        {
            ofile << A(i,j) << " "; 
        }
        
        ofile << std::endl;
        
    }
    
    // Close the output file
    ofile.close();

    return 0;
}

