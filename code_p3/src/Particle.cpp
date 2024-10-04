#include "utils.hpp"
#include "Particle.hpp"

#define EPSILON 1e-9

// constructor, the : instantiates member variables as opposed to local variables.
Particle::Particle(double m, double q, const arma::vec& r, const arma::vec& v) 
    : m(m), q(q), r(r), v(v) {
        // mostly empty constructor body, but check that the mass is not zero
        if (m < EPSILON) {
        std::cout << "Error: Mass of particle (" << m << ") cannot be near-zero or smaller." << std::endl;
        assert(m > EPSILON);
    }
} 

