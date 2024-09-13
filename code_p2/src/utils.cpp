// By including utils.hpp, we also include all 
// the headers included in utils.hpp
#include "utils.hpp"


// Calculating eigenvals and eigenvec analytical
td::pair<arma::vec, arma::mat> analytical(int N, double d, double a)
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


