#include "utils.hpp"
#include "Particle.hpp"

class PenningTrap {

    public:

        // Constructor
        PenningTrap(double B0_in, double V0_in, double d_in, bool interaction_in = true);

        // Attributes
        double B0;
        double V0;
        double d;
        std::vector<Particle> particles;
        bool interaction;

        // uninitialized
        double f;
        double w_v;

        // Methods
        void set_time_dependent_params(double f_in, double w_v_in);

        void add_particle(Particle p_in);

        void update_time_dependent_V0(double t);

        // External electric field at point r=(x,y,z)
        arma::vec external_E_field(arma::vec r);  

        // External magnetic field at point r=(x,y,z)
        arma::vec external_B_field(arma::vec r);  

        // Force on particle_i from particle_j
        arma::vec force_particle(int i, int j);

        // The total force on particle_i from the external fields
        arma::vec total_force_external(int i);

        // The total force on particle_i from the other particles
        arma::vec total_force_particles(int i);

        // The total force on particle_i from both external fields and other particles
        arma::vec total_force(int i);

        // Evolve the system one time step (dt) using Runge-Kutta 4th order
        void evolve_RK4(double dt);

        // Evolve the system one time step (dt) using Forward Euler
        void evolve_forward_Euler(double dt);

        void save_metadata(int num_particles, int num_timesteps, std::string filename);

        void save_to_file(const std::string& filename, double t);

};
// 