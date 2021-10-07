#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>
#include "example_header.hpp"

using namespace arma;
using namespace std;

class Particle
{
public:
    double charge;
    double mass;
    vec position;
    vec velocity;

    Particle(double q_in, double mass_in, vec position_in, vec velocity_in) //Constructor
    {
        charge = q_in;
        mass = mass_in;
        position = position_in;
        velocity = velocity_in;
    }

    void print_particle() // funtion
    {
        cout << "charge = " << charge << endl;
        cout << "mass = " << mass << endl;
        position.print("position = ");
        velocity.print("velocity = ");
    }
};


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
    void external_E_field(vec r)
    {

        vec external_E_field = vec(3, fill::zeros);

        external_E_field[0] = V0 * particles.r[0] / (d * d);

        external_E_field[1] = V0 * particles.r[1] / (d * d);

        external_E_field[2] = - V0 * 2 * particles.r[2] / (d * d);

    }

    // External magnetic field at point r=(x,y,z)
    void external_B_field(vec r)
    {
        vec external_B_field = vec(3, fill::zeros);

        external_B_field[2] = B0;

    }

    //Lorentz force
    void total_force_external(int i)
    {
      vec total_force_external = vec(3, fill::zeros);

      // Don't know if I use the right variable names !
      total_force_external[0] = p_in[i].q * external_E_field[0] + p_in[i].q * p_in[i].velocity_in * B0[2];

      total_force_external[1] = p_in[i].q * external_E_field[0] - p_in[i].q * p_in[i].velocity_in * B0[2];

      total_force_external[2] = 0;

    }

    // Force on particle_i from particle_j
    void total_force_particles(int i, int j)
    {
      double k_e = 1.38935333e5

      if (i != j){
        vec total_force_particles = k_e * q[i] / m[i]  //Got lost on how to implement

      }
    }


    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

};
