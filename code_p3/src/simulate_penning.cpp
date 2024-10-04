#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"

int simulate(std::vector<Particle> particles, PenningTrap trap, double dt, int timesteps, std::string filename, bool interactions = true) {

    for (const Particle& particle : particles) {
        trap.add_particle(particle);
    }

    if(!interactions) {
        trap.interaction = false; // might not work
    }

    std::ofstream ofile(filename, std::ios::trunc);
    ofile.close();

    trap.save_metadata(particles.size(), timesteps, filename);

    for (int i = 0; i < timesteps; i++) {
        trap.evolve_RK4(dt);
        trap.save_to_file(filename, i*dt);
    }

    return 0;
}


int main() {
    double timesteps = 51;
    double dt = 0.000001;
    arma::vec position = {20.0, 0.0, 20.0};
    arma::vec velocity = {0.0, 25.0, 0.0};

    PenningTrap trap(B0_converted, V0_converted, d_const);
    std::vector<Particle> particles;


    // Simulate the movement of a single particle in Penning trap for a total time of 50µs.
    // Make a plot of the motion in the direction as a function of time.
    particles.push_back(Particle(1.0, -1.0, {20.0, 0.0, 20.0}, {0.0, 25.0, 0.0})); // adding 1. electron
    simulate(particles, trap, dt, timesteps, "code_p3/data/one_particle_int.txt", true);

    // Simulate two particles in your Penning trap and make a plot of their motion in the xy-plane with and without particle interactions.
    particles.push_back(Particle(1.0, -1.0, {25.0, 25.0, 0.0}, {0.0, 40.0, 5.0})); // adding 2. electron
    simulate(particles, trap, dt, timesteps, "code_p3/data/two_particles_int.txt", true);

    // and without particle interactions
    simulate(particles, trap, dt, timesteps, "code_p3/data/two_particles_no_int.txt", false);

}