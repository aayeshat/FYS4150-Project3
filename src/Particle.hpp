#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

using namespace arma;

class Particle
{
public:
    double q; //charge
    double m; //mass
    vec r;    //position
    vec v;    //velocity

    Particle(double q_in, double m_in, vec r_in, vec v_in); //constructor
};

#endif