#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"


int main(){
    double B0 = 1;
    double V0 = 1;
    double d  = 1;
    double timesteps = 50;
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

    std::string filename = "data/testing.txt";
    std::ofstream ofile(filename, std::ios::trunc);
    ofile.close(); 

    for (int i = 0; i < timesteps; i++) {
        Trap.evolve_RK4(dt); 
        Trap.save_to_file(filename, i*dt, 2, timesteps); 
    }

    return 0;
}


//  g++ tests/test_penning.cpp src/utils.cpp src/Particle.cpp  src/PenningTrap.cpp -I include -o out -std=c++14 -I/opt/homebrew/Cellar/armadillo/14.0.2_1/include -L/opt/homebrew/Cellar/armadillo/14.0.2_1/lib -larmadillo && ./out

