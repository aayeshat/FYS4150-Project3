#include "Particle.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

class PenningTrap

{
public:
    double B0;
    double V0;
    double d;
    vector<Particle> particles;

    PenningTrap(double B0_in, double V0_in, double d_in) // Constructor
    {
        B0 = B0_in;
        V0 = V0_in;
        d = d_in;
    }
    void add_particle(Particle p_in)
    {
        particles.insert(particles.begin(), p_in);
    }

    // External electric field at point r=(x,y,z)
    vec external_E_field(vec r)
    {

        vec electricfield(3);
        return electricfield;
    }

    // External magnetic field at point r=(x,y,z)
    vec external_B_field(vec r)
    {
        vec magneticfield(3);
        return magneticfield ;
    }

    // Force on particle_i from particle_j
    vec force_particle(int i, int j)
    {
        vec force(2);

        return force;
    }
};

int main()
{
    PenningTrap pt(1., 3., 2);

    Particle p1(1., 40., vec(3, fill::ones), vec(3, fill::randu));
    pt.add_particle(p1);

    return 0;
}