#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"

void check_outside_region() {
    // task 9
    // Implement a check that sets the external E- and B-fields to zero in the region outside the trap. 
    // We’ll use the characteristic distance as a simple measure for the trap size, 
    // so what we need is a check that sets the external fields E and B to zero when |r| > d.
}


int main(){
    double B0 = B0_converted;
    double V0 = V0_converted;
    double d  = d_const;
    double tot_time = 0.00005;
    double timesteps = 30;
    double dt = tot_time/timesteps;
    
    double mass = 1.0;
    double charge = -1.0;
    arma::vec position = {1.0, 1.0, 1.0};
    arma::vec velocity = {1.0, 1.0, 1.0};

    double mass2 = 0.5;
    double charge2 = -2.0;
    arma::vec position2 = {0.0, 1.0, 1.0};
    arma::vec velocity2 = {1.0, 2.0, 1.0};

    bool interaction = false;
    int num_particles = 1;

    Particle electron1 = Particle(mass, charge, position, velocity);
    Particle electron2 = Particle(mass2, charge2, position2, velocity2);

    PenningTrap trap(B0, V0, d, interaction);
    trap.add_particle(electron1);
    //trap.add_particle(electron2);

    std::cout << "Bigger: " << B0*B0*1/40 << std::endl;
    std::cout << "Smaller: " << 4*V0/(d*d) << std::endl;


    //std::string filename = "data/testing_for_benjamin_one_particle.txt";
    //std::ofstream ofile(filename, std::ios::trunc);
    //ofile.close(); 

    //trap.save_metadata(num_particles, timesteps+1, filename);

    //for (int i = 0; i <= timesteps; i++) {
        //trap.save_to_file(filename, i*dt); 
        //trap.evolve_RK4(dt); 
    //}

    return 0;
}


//  g++ code_p3/src/simulate_penning.cpp code_p3/src/utils.cpp code_p3/src/Particle.cpp  code_p3/src/PenningTrap.cpp -I code_p3/include -o out -std=c++14 -I/opt/homebrew/Cellar/armadillo/14.0.2_1/include -L/opt/homebrew/Cellar/armadillo/14.0.2_1/lib -larmadillo && ./out

