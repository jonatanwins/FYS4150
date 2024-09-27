
#pragma once
#include "utils.hpp"

class Particle {
    
    public:
    
        // Constructor
        Particle(double m, double q, const arma::vec & r, const arma::vec& v);

        // Attributes
        double m;
        double q;
        arma::vec r;
        arma::vec v; 
};

