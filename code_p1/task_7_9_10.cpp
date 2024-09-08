#include "utils.hpp" 

int main() {

    int values[] = {10+1, 100+1, 1000+1, 10000+1, 100000+1, 1000000+1, 10000000+1};
    int num_n = sizeof(values) / sizeof(values[0]);
    std::vector<double> runtimes(num_n, 1.0); 
    std::vector<double> runtimes_special(num_n, 1.0); 
    
    for (int i = 0; i < num_n; ++i) {
        double n = values[i];

        double v_0 = 0.0;      
        double v_n_plus_1 = 0.0;

        std::vector<double> a(n, -1.0);
        std::vector<double> b(n+1, 2.0);
        std::vector<double> c(n, -1.0);

        std::vector<double> f(n+2);
        for (int i = 0; i <= n+1; i++) {
            f[i] = 100*exp(-10*(i/n));
        }

        // Exact solution
        std::vector<double> u = exactSolution(n);

        auto start = std::chrono::high_resolution_clock::now();

        // General method
        std::vector<double> v = solveTridiagonalMatrixEquation(a, b, c, f, v_0, v_n_plus_1);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> runtime = end - start;
        runtimes[i] = runtime.count();

        auto start_special = std::chrono::high_resolution_clock::now();

        // Special method
        std::vector<double> v_special = solveTridiagonalMatrixEquation_special(a, b, c, f, v_0, v_n_plus_1);   

        auto end_special = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> runtime_spec = end_special - start_special;
        runtimes_special[i] = runtime_spec.count();

        int int_n = static_cast<int>(n - 1);

        std::string filename = "data/v_and_u_n_" + std::to_string(int_n) + ".txt";
        w_file(filename, v, u);

    }

    std::vector<int> values_vector;
    for (int i = 0; i < num_n; ++i) {
        values_vector.push_back(values[i] - 1);
    }
    std::string filename_runtime = "data/runtime";
    w_file_one(filename_runtime, values_vector, runtimes, runtimes_special);

    return 0;

}



