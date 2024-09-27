#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"

PenningTrap(double B0_in, double V0_in, double d_in)
    : B0(B0_in), V0(V0_in), d(d_in), particles() {
}


// Add a particle to the trap
void add_particle(Particle p_in){

    particles.push_back(p_in);

} 

// External electric field at point r=(x,y,z)
arma::vec external_E_field(arma::vec r){

    double x = r.at(0);
    double y = r.at(1);
    double z = r.at(2);

    // E as the gradient of V
    double Ex = (V0*x) / (std::pow(d,2));
    double Ey = (V0*y) / (std::pow(d,2));
    double Ez = (-2*V0*z) / (std::pow(d,2));

    arma::vec E = arma::vec ({Ex, Ey, Ez});

    return E;
}

// External magnetic field at point r=(x,y,z)
arma::vec external_B_field(arma::vec r){

    arma::vec B = arma::vec ({0.0,0.0,this->B0});
    
    return B;
    
}  

// Force on particle_i from particle_j
arma::vec force_particle(int i, int j){
    // Columb force = k_e * q_1 * q_2 / norm(r)^3 * vector of distance between 

    arma::vec r_i = particles.at(i).r;
    arma::vec r_j = particles.at(j).r;

    double q_i = particles.at(i).q;
    double q_i = particles.at(j).q;

    arma::vec r_diff = r_i - r_j;

    double r_diff_norm = arma::norm(r_diff);

    arma::vec force_ij = (k_e * q_i * q_j)/std::pow(r_diff_norm,3) * r_diff;

    return force_ij;
    
}

// The total force on particle_i from the external fields
arma::vec total_force_external(int i){
    // F = qE + qv x B

    double q_i = particles.at(i).q;
    double v_i = particles.at(i).v;
    double r_i = particles.at(i).r;

    arma::vec E_field = external_E_field(r_i);
    arma::vec B_field = external_B_field(r_i);

    arma::vec external_forces = q_i * E_field + arma::cross(q_i * v_i, B_field);

    return external_forces;
    
}

// The total force on particle_i from the other particles
arma::vec total_force_particles(int i){

    arma::vec force_particles = arma::vec({0.0, 0.0, 0.0});
    // loop over all other particles j
    for (j = 0, j = particles.size(), j++){

        arma::double force_ij = force_particle(i,j);
        // add to force_particles
        force_particles += force_ij;
    }

    return force_particles;
    
}

// The total force on particle_i from both external fields and other particles
arma::vec total_force(int i){

    arma::vec force_particles = total_force_particles(i);
    arma::vec external_forces = total_force_external(i);

    arma::vec total_forces = total_force_particles + total_force_external;

    return total_forces;
    
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void evolve_RK4(double dt){
    // temporary copy of all the particles in the Penning trap
    std::vector particles_copy = particles;

    // init 
    arma::vec k_r_1, k_r_2, k_r_3, k_r_4;
    arma::vec k_v_1, k_v_2, k_v_3, k_v_4;

    // k1 
    for (j = 0, j = particles.size(), j++){
        arma::vec k_r_1 = dt * particles_copy.at(j).v;
        arma::vec k_v_1 = dt * total_forces(j)/particles_copy.at(j).m;
    } ;

    // update after k1 
    for (size_t j = 0; j < particles.size(); j++) {
        particles.at(j).r = particles_copy.at(j).r + 0.5 * dt * k_r_1;
        particles.at(j).v = particles_copy.at(j).v + 0.5 * dt * k_v_1;
    } ;

    // need to update time??? t += dt/2

    // k2 
    for (j = 0, j = particles.size(), j++){
        arma::vec k_r_2 = dt * particles.at(j).v;
        arma::vec k_v_2 = dt * total_forces(j)/particles.at(j).m;
    } ;

    // update after k2
    for (size_t j = 0; j < particles.size(); j++) {
        particles.at(j).r = particles_copy.at(j).r + 0.5 * dt * k_r_2;
        particles.at(j).v = particles_copy.at(j).v + 0.5 * dt * k_v_2;
    } ;
    
    // need to update time??? t += dt/2

    // k3 
    for (j = 0, j = particles.size(), j++){
        arma::vec k_r_3 = dt * particles.at(j).v;
        arma::vec k_v_3 = dt * total_forces(j)/particles.at(j).m;
    } ;

    // update after k3
    for (size_t j = 0; j < particles.size(); j++) {
        particles.at(j).r = particles_copy.at(j).r + dt * k_r_3;
        particles.at(j).v = particles_copy.at(j).v + dt * k_v_3;
    }  ;

    // need to update time??? t += dt

    // k4
    for (j = 0, j = particles.size(), j++){
        arma::vec k_r_4 = dt * particles.at(j).v;
        arma::vec k_v_4 = dt * total_forces(j)/particles.at(j).m;
    } ; 

    // final update 
    for (j = 0, j = particles.size(), j++){
        particles.at(j).r = particles_copy.at(j).r + (1/6) * (k_r_1 + 2*k_r_2 + 2*k_r_3 + k_r_4);
        particles.at(j).v = particles_copy.at(j).v + (1/6) * (k_v_1 + 2*k_v_2 + 2*k_v_3 + k_v_4);
    } ; 
    
} 

// Evolve the system one time step (dt) using Forward Euler
void evolve_forward_Euler(double dt){

    for (j = 0, j = particles.size(), j++){
        particles.at(j).r = particles.at(j).r + dt * particles.at(j).v;
        particles.at(j).v = particles.at(j).v + dt * total_forces(j)/particles.at(j).m;
    } ; 

}  
        