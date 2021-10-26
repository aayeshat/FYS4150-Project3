#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp_

#include "../Particle.hpp"
#include <armadillo>

using namespace arma;
using namespace std;

class PenningTrap
{
public:
    double B0;
    double V0;
    double d;
    int n;

    double ke = 1.38935333e5;
    vector<Particle> particles;

    mat r_step;
    mat v_step;
    bool interaction = true;

    PenningTrap(double B0_in, double V0_in, double d_in, int n_in);

    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    vec external_E_field(int i);

    // External magnetic field at point r=(x,y,z)
    vec external_B_field(int i);

    // coulombic interactions Force on particle_i from particle_j
    vec force_particle(int i, int j);

    //total force on particle i due to external fields

    vec total_force_external(int i);

    //force on particle i from other particles
    vec total_force_particles(int i);

    //force on particle due to fields and particles
    vec total_force(int i);

    void evolve_RK4(double dt);

    void evolve_forward_Euler(double dt);
};

#endif