#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"
#include <vector>
#include <iostream>

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool interaction_in)
    : B0(B0_in), V0(V0_in), d(d_in), particles(), interaction(interaction_in) {
}

void PenningTrap::set_time_dependent_params(double f_in, double w_v_in) {
    f = f_in;
    w_v = w_v_in;
}

void PenningTrap::update_time_dependent_V0(double t) {
    double V0_dep = V0_converted*(1 + f*std::cos(w_v * t));
    V0 = V0_dep;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in){

    particles.push_back(p_in);

} 

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){

    double x = r.at(0);
    double y = r.at(1);
    double z = r.at(2);

    // E as the gradient of V
    double Ex = (this->V0*x) / (std::pow(this->d,2));
    double Ey = (this->V0*y) / (std::pow(this->d,2));
    double Ez = (-2*this->V0*z) / (std::pow(this->d,2));

    arma::vec E = arma::vec ({Ex, Ey, Ez});

    return E;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r){

    arma::vec B = arma::vec ({0.0,0.0,this->B0});
    
    return B;
    
}  

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j){
    // Columb force = k_e * q_1 * q_2 / norm(r)^3 * vector of distance between 

    arma::vec force_ij = {0.0, 0.0, 0.0};

    if (i != j) {
        arma::vec r_i = particles.at(i).r;
        arma::vec r_j = particles.at(j).r;

        double q_i = particles.at(i).q;
        double q_j = particles.at(j).q;

        arma::vec r_diff = r_i - r_j;

        double r_diff_norm = arma::norm(r_diff);

        force_ij = (k_e * q_i * q_j)/std::pow(r_diff_norm,3) * r_diff;
    }
        
    return force_ij;
    
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
    // F = qE + qv x B

    double q_i = particles.at(i).q;
    arma::vec v_i = particles.at(i).v;
    arma::vec r_i = particles.at(i).r;

    arma::vec E_field = PenningTrap::external_E_field(r_i);
    arma::vec B_field = PenningTrap::external_B_field(r_i);

    arma::vec external_forces = q_i * E_field + arma::cross(q_i * v_i, B_field);

    return external_forces;
    
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){

    arma::vec force_particles = arma::vec({0.0, 0.0, 0.0});

    // loop over all other particles j
    for (int j = 0; j < particles.size(); j++){

        arma::vec force_ij = PenningTrap::force_particle(i,j);

        // add to force_particles
        force_particles += force_ij;
    }

    return force_particles;
    
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i){

    arma::vec external_forces = PenningTrap::total_force_external(i);
    
    if (arma::norm(particles.at(i).r) > d_const) {
        external_forces = {0.0, 0.0, 0.0};
    }

    arma::vec force_particles = {0.0, 0.0, 0.0};
    if (this->interaction == true) {
        //std::cout << "interaction is true" << std::endl;
        force_particles = PenningTrap::total_force_particles(i);
    }

    arma::vec total_forces =  external_forces + force_particles;

    return total_forces;
    
}


void PenningTrap::evolve_RK4(double dt){
    // Define the number of dimensions and particles
    int particle_dim = 3;  // assuming 3D space (x, y, z)
    int num_particles = particles.size();  // total number of particles

    // Initialize k_r and k_v for each particle
    std::vector<arma::Col<double>> k_r_1(num_particles), k_r_2(num_particles), k_r_3(num_particles), k_r_4(num_particles);
    std::vector<arma::Col<double>> k_v_1(num_particles), k_v_2(num_particles), k_v_3(num_particles), k_v_4(num_particles);

    // Copy the particles
    std::vector<Particle> particles_copy = particles;

    // k1
    for (int j = 0; j < num_particles; j++){
        k_r_1[j] = dt * particles_copy[j].v;  // particles_copy[j].v is arma::Col<double>
        k_v_1[j] = dt * PenningTrap::total_force(j) / particles_copy[j].m;  // total_force returns arma::Col<double>
    }

    // Update after k1
    for (int j = 0; j < num_particles; j++) {
        particles[j].r += 0.5 * k_r_1[j];  // r is arma::Col<double>
        particles[j].v += 0.5 * k_v_1[j];  // v is arma::Col<double>
    }

    // k2
    for (int j = 0; j < num_particles; j++){
        k_r_2[j] = dt * particles[j].v;
        k_v_2[j] = dt * PenningTrap::total_force(j) / particles[j].m;
    }

    // Update after k2
    for (int j = 0; j < num_particles; j++) {
        particles[j].r += 0.5 * k_r_2[j];
        particles[j].v += 0.5 * k_v_2[j];
    }

    // k3
    for (int j = 0; j < num_particles; j++){
        k_r_3[j] = dt * particles[j].v;
        k_v_3[j] = dt * PenningTrap::total_force(j) / particles[j].m;
    }

    // Update after k3
    for (int j = 0; j < num_particles; j++) {
        particles[j].r += k_r_3[j];
        particles[j].v += k_v_3[j];
    }

    // k4
    for (int j = 0; j < num_particles; j++){
        k_r_4[j] = dt * particles[j].v;
        k_v_4[j] = dt * PenningTrap::total_force(j) / particles[j].m;
    }

    // Final update
    for (int j = 0; j < num_particles; j++){
        particles[j].r = particles_copy[j].r + (1.0/6.0) * (k_r_1[j] + 2.0*k_r_2[j] + 2.0*k_r_3[j] + k_r_4[j]);
        particles[j].v = particles_copy[j].v + (1.0/6.0) * (k_v_1[j] + 2.0*k_v_2[j] + 2.0*k_v_3[j] + k_v_4[j]);
    }
}


// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){

    for (int j = 0; j < particles.size(); j++){
        particles.at(j).r = particles.at(j).r + dt * particles.at(j).v;
        particles.at(j).v = particles.at(j).v + dt * PenningTrap::total_force(j)/particles.at(j).m;
    } 

}  

void PenningTrap::save_metadata(int num_particles, int num_timesteps, std::string filename) {
    std::ofstream ofile(filename, std::ios::app);

    if (!ofile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return; 
    }

    if (ofile.tellp() == 0) {
        ofile << num_particles << " " << num_timesteps << std::endl; 
    }
    else {
        throw std::runtime_error("Error: Cannot add metadata to nonempty file, " + filename + " is not empty.");
    }
}

void PenningTrap::save_to_file(const std::string& filename, double t) {
    
    // open file in append mode
    std::ofstream ofile(filename, std::ios::app); 
    //  TODO clean up what chatgpt has written 

    // print error message if problems with opening file
    if (!ofile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return; 
    }

    // write time to file
    ofile << t << std::endl;

    // write particle data to file
    int precision = 10; 
    for (int i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];

        for (int k = 0; k < particle.r.n_elem; k++) {
            ofile << std::setprecision(precision) << std::scientific << particle.r(k) << " "; 
        }

        for (int k = 0; k < particle.v.n_elem; k++) {
            ofile << std::setprecision(precision) << std::scientific << particle.v(k) << " "; 
        }

        ofile << std::endl; 
    }

    ofile.close(); 
}
