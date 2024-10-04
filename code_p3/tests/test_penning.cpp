#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"



// -------------------------------------- Constants ------------------------------------
// Coulomb constant 
// double k_e = 1.38935333e5; // (μm)^3 / (μs)^2 e^2

// Magnetic field strength (Tesla) and electric potential (Volt)
// double T_unit = 9.64852558e1; // μm / (μs * e)
// double V_unit = 9.64852558e7; // μm^2 / (μs^2 * e)

// Penning trap configuration
// double B0 = 1.00; // Tesla
// double B0_converted = B0 * T_unit; // μm / (μs * e)

// double V0 = 25.0e-3; // 25.0 mV = 25.0e-3 V
// double V0_converted = V0 * V_unit; // μm^2 / (μs^2 * e)

// double d = 500.0; // μm

// Ratio of V0/d^2
// double V0_d2 = V0_converted / std::pow(d, 2);

// -------------------------------------------------------------------------------------



int main(){
    double B0 = 1;
    double V0 = 1;
    double d  = 1;
    double timesteps = 2;
    double dt = 0.01;
    
    double mass = 1.0;
    double charge = -1.0;
    arma::vec position = {1.0, 1.0, 1.0};
    arma::vec velocity = {1.0, 1.0, 1.0};

    double mass2 = 0.5;
    double charge2 = -2.0;
    arma::vec position2 = {0.0, 1.0, 1.0};
    arma::vec velocity2 = {1.0, 2.0, 1.0};

    Particle electron1 = Particle(mass, charge, position, velocity);
    Particle electron2 = Particle(mass2, charge2, position2, velocity2);

    PenningTrap Trap(B0, V0, d);
    Trap.add_particle(electron1);
    Trap.add_particle(electron2);

    for (int i = 0; i < timesteps; i++) {
        Trap.evolve_forward_Euler(dt);
    }

    return 0;
}


//  g++ tests/test_penning.cpp src/utils.cpp src/Particle.cpp  src/PenningTrap.cpp -I include -o out -std=c++14 -I/opt/homebrew/Cellar/armadillo/14.0.2_1/include -L/opt/homebrew/Cellar/armadillo/14.0.2_1/lib -larmadillo && ./out

