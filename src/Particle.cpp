#include "Particle.hpp"

#include <armadillo>

using namespace arma;

Particle::Particle(double q_in, double m_in, vec r_in, vec v_in) //constructor
{
    q = q_in;
    m = m_in;
    r = r_in;
    v = v_in;
}