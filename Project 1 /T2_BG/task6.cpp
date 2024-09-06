#include "utils.hpp" // Assuming this file contains necessary declarations

// Forward substitution function

std::pair<std::vector<double>, std::vector<double>> forwardSubstitution(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1 ) {

    int n = b.size();
    double h = 1./n;

    std::vector<double> g(n+2);

    g[1] = h*h*f[1] + v_0;
    g[n] = h*h*f[n] + v_n_plus_1;

    for (int i = 2; i <= n-1; i++) {
        g[i] = h*h*f[i];
    }

    std::vector<double> btilde(n+2);
    std::vector<double> gtilde(n+2);

    btilde[1] = b[1];
    gtilde[1] = g[1];
 
    for (int i = 2; i <= n; i++) {
        btilde[i] = b[i] - (a[i]/btilde[i-1])*c[i-1];
        gtilde[i] = g[i] - (a[i]/btilde[i-1])*gtilde[i-1];
    }

    return {gtilde, btilde};   
}

// Backward substitution function

std::vector<double> backwardSubstitution(const std::vector<double>& gtilde, const std::vector<double>& btilde, const std::vector<double>& c, double v_0, double v_n_plus_1) {

    int n = c.size()+1;
    std::vector<double> v(n+2);

    v[0] = v_0;
    v[n+1] = v_n_plus_1;

    v[n] = gtilde[n]/btilde[n];

    for (int i = n-1; i >= 1; i--) {
        v[i] = (gtilde[i] - v[i+1]*c[i])/btilde[i];
    }

    return v;
}

// Function to solve the tridiagonal matrix equation

std::vector<double> solveTridiagonalMatrixEquation(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& f, double v_0, double v_n_plus_1) {

    std::pair<std::vector<double>, std::vector<double>> forw = forwardSubstitution(a, b, c, f, v_0, v_n_plus_1);
    std::vector<double> gtilde = forw.first;
    std::vector<double> btilde = forw.second;
    return backwardSubstitution(gtilde, btilde, c, v_0, v_n_plus_1);
}

// Exact solution : u(x) = 1 - (1-e^{-10})x - e^{-10x}

std::vector<double> exactSolution(double n) {
    
    std::vector<double> u(n+2);
    for (int i = 0; i <= n+1; i++) {
        u[i] = 1 - (1 - exp(-10))*(i/n) - exp(-10*(i/n));
    }
    return u;
}

int main() {
    int n = 100;           
    double v_0 = 0.0;      
    double v_n_plus_1 = 0.0;  

    std::vector<double> a(n-1, -1.0);
    std::vector<double> b(n, 2.0);
    std::vector<double> c(n-1, -1.0);

    std::vector<double> f(n+2);
    for (int i = 1; i <= n; i++) {
        f[i] = 100*exp(-10*(i/n));
    }

    f[0] = v_0;
    f[n+1] = v_n_plus_1;

    std::vector<double> v = solveTridiagonalMatrixEquation(a, b, c, f, v_0, v_n_plus_1);
    std::vector<double> u = exactSolution(n);

    // print v versus u 

    for (int i = 0; i <= n+2; i++) {
        std::cout << v[i] << " " << u[i] << std::endl;
    }
    return 0;

}

