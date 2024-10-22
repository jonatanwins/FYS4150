#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"

void create_directories() {

    std::vector<std::string> dirs = {"data", "plots"};
    
    for (const std::string& dir : dirs) {
        std::filesystem::create_directories(dir);
    }
}

void simulate(std::vector<Particle> particles, PenningTrap trap, double dt, int timesteps, 
        std::string filename, bool interactions = true, std::string method = "RK4", bool time_dependent_E = false, int sample_rate = 1) {

    if (sample_rate == -1) {
        sample_rate = timesteps;
    }

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

        if (i % sample_rate == 0) {
            trap.save_to_file(filename, i*dt);
        }
        
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

    double time = 50;
    double dt = 0.01;
    int timesteps = time/dt;


    // Simulate the movement of a single particle in Penning trap for a total time of 50µs.
    simulate(particles, trap, dt, timesteps, "data/one_particle_int_RK4.txt", true, "RK4", false);
    simulate(particles, trap, dt, timesteps, "data/one_particle_int_FE.txt", true, "FE", false);

    // Simulate two particles in your Penning trap and make a plot of their motion in the xy-plane with and without particle interactions.
    particles.push_back(Particle(40.078, 1.0, {25.0, 25.0, 0.0}, {0.0, 40.0, 5.0})); // adding 2. proton
    simulate(particles, trap, dt, timesteps, "data/two_particles_int_RK4.txt",true,"RK4", false);
    simulate(particles, trap, dt, timesteps, "data/two_particles_int_FE.txt",true,"FE", false);

    // and without particle interactions
    simulate(particles, trap, dt, timesteps, "data/two_particles_no_int_RK4.txt", false, "RK4", false);
    simulate(particles, trap, dt, timesteps, "data/two_particles_no_int_FE.txt", false, "FE", false);

    //  1 particle and simulation time 50µs. Run the simulation four times, using 4000, 8000, 16000, 32000 timesteps
    // Do the same using the forward Euler method.

    particles.pop_back(); // remove second particle
    for (const int& n : {4000, 8000, 16000, 32000}) {
        std::ostringstream filename;
        filename << "data/one_particle_no_int_n=" << n << "_RK4.txt";
        simulate(particles, trap, time/n, n, filename.str(), false, "RK4", false);
        
        filename.str("");
        
        filename << "data/one_particle_no_int_n=" << n << "_FE.txt";
        simulate(particles, trap, time/n, n, filename.str(), false, "FE", false);
    
    }
}

void simulate_traps_time_dependent_E(std::vector<Particle> particles, PenningTrap trap, bool interactions, int sample_rate) {

    double time = 500;
    double dt = 0.01;
    int timesteps = time/dt;

    int size_particles = particles.size();
    std::ostringstream filename;

    if (interactions && sample_rate == 1000) {
        filename << "data/int_" << "f_" << trap.f << "_w_v_" << trap.w_v << "_.txt";
    }
    else if (interactions && sample_rate == 10) {
        filename << "data/int_" << "f_" << trap.f << "_w_v_" << trap.w_v << "_full.txt";
    }
    else {
        filename << "data/no_int_" << "f_" << trap.f << "_w_v_" << trap.w_v << "_.txt";
    }
    

    simulate(particles, trap, dt, timesteps, filename.str(), interactions, "RK4", true, sample_rate);
    std::cout << "Completed simulation " << filename.str() << std::endl;
}

void simulate_arbitrary_particles(std::vector<Particle> particles, PenningTrap trap, int number_of_particles, int sample_rate, bool interactions = false,
                                    double f_start = 0.1, double f_stop = 0.7, double f_step = 0.3, 
                                    double w_v_start = 2.18, double w_v_stop = 2.32, double w_v_step = 0.02
                                    ) {
    arma::arma_rng::set_seed(4150); // set seed for reproducability, FYS4150
    
    for (int i = 0; i < number_of_particles; i++) {
        particles.push_back(Particle(40.078, 1.0, arma::vec(3).randn() * 0.1 * d_const, arma::vec(3).randn()*0.1*d_const)); // adding more protons
    }


    if (f_step != 0.0 && w_v_step != 0.0) {
        
    for (double f = f_start; f <= f_stop; f += f_step) {
        for (double w_v = w_v_start; w_v <= w_v_stop; w_v += w_v_step) { 
            std::cout << "f: " << f << " w_v: " << w_v << std::endl;
            trap.set_time_dependent_params(f, w_v);
            simulate_traps_time_dependent_E(particles, trap, interactions, sample_rate);
            }
        }
    }

    else if (f_step == 0.0 && w_v_step != 0.0) {
        for (double w_v = w_v_start; w_v <= w_v_stop; w_v += w_v_step) { 
            std::cout << "w_v: " << w_v << std::endl;
            trap.set_time_dependent_params(0.0, w_v);
            simulate_traps_time_dependent_E(particles, trap, interactions, sample_rate);
        }
    }

    else if (f_step != 0.0 && w_v_step == 0.0) {
        for (double f = f_start; f <= f_stop; f += f_step) {
            std::cout << "f: " << f << std::endl;
            trap.set_time_dependent_params(f, 0.0);
            simulate_traps_time_dependent_E(particles, trap, interactions, sample_rate);
        }
    }

    else {
        trap.set_time_dependent_params(f_start, w_v_start);
        simulate_traps_time_dependent_E(particles, trap, interactions, sample_rate);
        }
}
                                    


int main() {

    // initialization 
    create_directories();
    PenningTrap trap(B0_converted, V0_converted, d_const);
    std::vector<Particle> particles;

    // --- for the simulation of one and two particles with constant E, with and without interactions ---
    simulate_traps_constant_E(particles, trap);

    // --- for the grid search (no interactions) ---> expected file size = 11 kB 
    //simulate_arbitrary_particles(particles, trap, 100, 50001, false, 0.1, 0.7, 0.3, 0.2, 2.52, 0.02);

    // --- for the exploration of w_v found from the grid search (with interactions) ---> expected file size = 541 kB 
    //simulate_arbitrary_particles(particles, trap, 100, 1000, true, 0.4, 0.4, 0.0, 0.67, 0.67, 0.00);
    //simulate_arbitrary_particles(particles, trap, 100, 1000, true, 0.4, 0.4, 0.0, 1.36, 1.4, 0.01);
    //simulate_arbitrary_particles(particles, trap, 100, 1000, true, 0.4, 0.4, 0.0, 2.2, 2.24, 0.01);

    // --- smaller sample rate for the best w_v (with interactions) ---> expected file size = 52.9 MB
    //simulate_arbitrary_particles(particles, trap, 100, 10, true, 0.4, 0.4, 0.0, 0.68, 0.68, 0.01);
    //simulate_arbitrary_particles(particles, trap, 100, 10, true, 0.4, 0.4, 0.0, 1.38, 1.38, 0.01);
    //simulate_arbitrary_particles(particles, trap, 100, 10, true, 0.4, 0.4, 0.0, 2.21, 2.21, 0.01);

    
}

