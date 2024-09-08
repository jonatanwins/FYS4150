#include "utils.hpp"


int main()
{
    std::vector<double> x(100); // vector x lenght 100
    std::vector<double> u(x.size());

    double step = 1.0 / (x.size() -1);
    for (size_t i = 0; i < x.size(); i++)
    {
        x[i] = i * step;
        u[i] = 1 - (1 - exp(-10)) * x[i] - exp(-10 * x[i]);
    }


    // Write data to output file
    std::string filename = "data/output.txt";
    w_file(filename, x, u);
    
    return 0;
}

