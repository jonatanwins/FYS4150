#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"

int main(){
    double B0 = 1;
    double V0 = 1;
    double d  = 1; 
    PenningTrap Trap(B0, V0, d);
    // std::cout << Trap.B0 << std::endl;
    

    return 0;
}


//  g++ tests/test_penning.cpp src/utils.cpp src/Particle.cpp  src/PenningTrap.cpp -I include -o out -std=c++14 -I/opt/homebrew/Cellar/armadillo/14.0.2_1/include -L/opt/homebrew/Cellar/armadillo/14.0.2_1/lib -larmadillo && ./out

