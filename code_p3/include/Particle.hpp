
#pragma once
#include "utils.hpp"

class Particle {
    
    public:
    
        // Constructor
        Particle(double m, double q, const arma::vec & pos, const arma::vec& vel);

        // Attributes
        double mass;
        double charge;
        arma::vec position;
        arma::vec velocity; 
};

