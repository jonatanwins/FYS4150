#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"

int simulate(str filename, bool interactions, int no_particles) {
    
    Trap.add_particle(electron);

    for (int i = 0; i < timesteps; i++) {
        Trap.evolve_RK4(dt);
        Trap.save_to_file(filename, i*dt);
    }

    return 0;
}


int main() {
    arma::vec position = {20.0, 0.0, 20.0};
    arma::vec velocity = {0.0, 25.0, 0.0};
    Particle electron = Particle(1.0, -1.0, position, velocity);

    double B0 = 1.0;
    double V0 = 1.0;
    double d  = 1.0;
    double timesteps = 50;
    double dt = 0.01;
    PenningTrap Trap(B0, V0, d);

    
    simulate("single_particle.txt", );
    // simulate_two_particles();
}