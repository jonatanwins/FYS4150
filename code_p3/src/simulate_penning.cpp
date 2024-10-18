#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"
#include <filesystem>
#include <armadillo>

void create_directories() {
    std::vector<std::string> dirs = {"data", "plots"};
    
    for (const std::string& dir : dirs) {
        std::filesystem::create_directories(dir);
    }
}

void simulate(std::vector<Particle> particles, PenningTrap trap, double dt, int timesteps, 
        std::string filename, bool interactions = true, std::string method = "RK4", bool time_dependent_E = false) {

    for (const Particle& particle : particles) {
        trap.add_particle(particle);
    }

    if(!interactions) {
        trap.interaction = false; // might not work
    }

    std::ofstream ofile(filename, std::ios::trunc);
    ofile.close();

    trap.save_metadata(particles.size(), timesteps+1, filename);

    if (method == "RK4") {
        for (int i = 0; i < timesteps; i++) {
            if (time_dependent_E) {
                trap.update_time_dependent_V0(i*dt);
            }
        trap.save_to_file(filename, i*dt);
        trap.evolve_RK4(dt);
        }
        trap.save_to_file(filename, timesteps*dt);
    }
    
    else if(method == "FE") {
        for (int i = 0; i < timesteps; i++) {
        trap.save_to_file(filename, i*dt);
        trap.evolve_forward_Euler(dt);
        }
        trap.save_to_file(filename, timesteps*dt);
    }

    else {
        std::cout << "Error: method for numerical differentiation must be RK4 or FE." << std::endl; 
    }
}

void simulate_traps_constant_E(std::vector<Particle> particles, PenningTrap trap) {

    particles.push_back(Particle(40.078, 1.0, {20.0, 0.0, 20.0}, {0.0, 25.0, 0.0})); // adding 1. proton
    particles.push_back(Particle(40.078, 1.0, {25.0, 25.0, 0.0}, {0.0, 40.0, 5.0})); // adding 2. proton

    double time = 500;
    double dt = 0.01;
    int timesteps = time/dt;

    // Simulate the movement of a single particle in Penning trap for a total time of 50µs.
    // Make a plot of the motion in the direction as a function of time.
    simulate(particles, trap, dt, timesteps, "data/one_particle_int_RK4.txt", true);

    // Simulate two particles in your Penning trap and make a plot of their motion in the xy-plane with and without particle interactions.
    particles.push_back(Particle(40.078, 1.0, {25.0, 25.0, 0.0}, {0.0, 40.0, 5.0})); // adding 2. proton
    simulate(particles, trap, dt, timesteps, "data/two_particles_int_RK4.txt", true);

    // and without particle interactions
    simulate(particles, trap, dt, timesteps, "data/two_particles_no_int_RK4.txt", false);

    //  1 particle and simulation time 50µs. Run the simulation four times, using 4000, 8000, 16000, 32000 timesteps
    // Do the same using the forward Euler method.
    particles.pop_back(); // remove second particle
    for (const int& n : {4000, 8000, 16000, 32000}) {
        std::ostringstream filename;
        filename << "data/one_particle_no_int_n=" << n << "_RK4.txt";
        simulate(particles, trap, timesteps/n, n, filename.str(), true, "RK4");
        
        filename.str("");
        
        filename << "data/one_particle_no_int_n=" << n << "_FE.txt";
        simulate(particles, trap, timesteps/n, n, filename.str(), true, "FE");
    
    }
}

void simulate_traps_time_dependent_E(std::vector<Particle> particles, PenningTrap trap, bool interactions) {
    double time = 500;
    double dt = 0.01;
    int timesteps = time/dt;

    int size_particles = particles.size();
    std::ostringstream filename;

    if (interactions) {
        filename << "data/int_" << "f_" << trap.f << "_w_v_" << trap.w_v << "_.txt";
    }
    else {
        filename << "data/no_int_" << "f_" << trap.f << "_w_v_" << trap.w_v << "_.txt";
    }
    

    simulate(particles, trap, dt, timesteps, filename.str(), interactions, "RK4", true);
    std::cout << "Completed simulation " << filename.str() << std::endl;
}

void simulate_arbitrary_particles(std::vector<Particle> particles, PenningTrap trap, int number_of_particles, bool interactions = false,
                                    double f_start = 0.1, double f_stop = 0.7, double f_step = 0.3, 
                                    double w_v_start = 2.18, double w_v_stop = 2.32, double w_v_step = 0.2
                                    ) {
    arma::arma_rng::set_seed(4150); // set seed for reproducability, FYS4150
    
    for (int i = 0; i < number_of_particles; i++) {
        particles.push_back(Particle(40.078, 1.0, arma::vec(3).randn() * 0.1 * d_const, arma::vec(3).randn()*0.1*d_const)); // adding more protons
    }

    for (double f = f_start; f <= f_stop; f += f_step) {
        for (double w_v = w_v_start; w_v <= w_v_stop; w_v += w_v_step) { 
            trap.set_time_dependent_params(f, w_v);
            simulate_traps_time_dependent_E(particles, trap, interactions);
        }
    }
}
                                    


int main() {
    // intialization
    create_directories();
    PenningTrap trap(B0_converted, V0_converted, d_const);
    std::vector<Particle> particles;

    // choose simulation
    simulate_arbitrary_particles(particles, trap, 100, true, 
    0.4, 0.4, 0.01, 
    2.15, 2.25, 0.01);
    
}
    
    
