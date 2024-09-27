#include "utils.hpp"

// -------------------------------------- Constants ------------------------------------
// Coulomb constant 
double ke = 1.38935333e5; // (μm)^3 / (μs)^2 e^2

// Magnetic field strength (Tesla) and electric potential (Volt)
double T_unit = 9.64852558e1; // μm / (μs * e)
double V_unit = 9.64852558e7; // μm^2 / (μs^2 * e)

// Penning trap configuration
double B0 = 1.00; // Tesla
double B0_converted = B0 * T_unit; // μm / (μs * e)

double V0 = 25.0e-3; // 25.0 mV = 25.0e-3 V
double V0_converted = V0 * V_unit; // μm^2 / (μs^2 * e)

double d = 500.0; // μm

// Ratio of V0/d^2
double V0_d2 = V0_converted / std::pow(d, 2);

// -------------------------------------------------------------------------------------


int main(){
      
    
    return 0;
}

