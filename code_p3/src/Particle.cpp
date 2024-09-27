#include "utils.hpp"
#include "Particle.hpp"

// constructor, the : instantiates member variables as opposed to local variables.
Particle::Particle(double m, double q, const arma::vec& r, const arma::vec& v) 
    : m(m), q(q), r(r), v(v) {
        // empty constructor body as the particles are only instantiated
} 

