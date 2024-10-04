#include "utils.hpp"
#include "PenningTrap.hpp"
#include "Particle.hpp"

void init_particle() {
    double mass = 1.0;
    double charge = -1.0;
    arma::vec position = {0.0, 0.0, 0.0};
    arma::vec velocity = {1.0, 0.0, 0.0};

    Particle electron = Particle(mass, charge, position, velocity);
    Particle antimatter = Particle(-7.9, -1.0, position, velocity);

    std::cout << "test for particle initialization passed" << std::endl;  
}

int main() {
    init_particle();
} 