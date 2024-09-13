#include "utils.hpp" 



int main()
{
    // Define a 4x4 matrix A
    arma::mat A_test = {
        {1, 0, 0, 0.5},
        {0, 1, -0.7, 0},
        {0, -0.7, 1, 0},
        {0.5, 0, 0, 1}
    };
}

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l)
{
    assert(A.is_square());
    assert(A.is_symmetric());

    double max_value = 0;

    /* Trenger bare sjekke øvre halvdel */
    for (size_t i = 0; i < A.n_rows; i++) {
        for (size_t j=i+1; j < A.n_cols; j++) {
            if (std::abs(A(i,j)) > max_value) {
                max_value = A(i,j);
                k = i;
                l = j;
            }
        }
    }
    return max_value;
    
    // for (size_t i = 0; i < A.n_rows; i++)
    // {
    //     for (size_t j = 1; i < A.n_cols; i++)
    //     {
    //         if (A(i,j) > A(i,j-1))
    //         {
    //             k = i
    //             l = j
                
    //         }
            
    //     }
        
    // }
    
}