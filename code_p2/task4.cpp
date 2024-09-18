#include "utils.hpp" 

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


int main() {

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

    double eps = 1e-4;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1e5;
    int iterations = 0;
    bool converged = true;

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
    // std::cout << "Eigenvalues: "<< eigenvalues << std::endl;
    // std::cout << "Eigenvectors: "<< eigenvectors << std::endl;
    // std::cout << "Iterations: "<< iterations << std::endl;
    // std::cout << "Converged: "<< converged << std::endl;

    eigenvalues.print("Jacobi Rotation Eigvalues");
    eigenvectors.print("Jacobi Rotation Eigenvectors");
}
    