#include "utils.hpp" 

int main()
{
    
    arma::mat A_test = {
        {1, 0, 0, 0.5},
        {0, 1, -0.7, 0},
        {0, -0.7, 1, 0},
        {0.5, 0, 0, 1}
    };

    int k = 0; 
    int l = 0;

    // Printing out largest off-diagonal values and indices with max_offdiag_symmetric(), see utils.cpp in src/
    std::cout << "max diagonal offset:" << max_offdiag_symmetric(A_test, k, l) << " k: " << k << " l: " << l << std::endl;
    return 0;
}