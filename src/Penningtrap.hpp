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
    int n;
    double ke = 1.38935333e5;
    vector<Particle> particles;

    PenningTrap(double B0_in, double V0_in, double d_in, int n_in) // Constructor
    {
        B0 = B0_in; //magnetic filed
        V0 = V0_in; //Applied Potential
        d = d_in;   //dimension
        n = n_in;
    }
    void add_particle(Particle p_in)
    {
        //particles.insert(particles.begin(), p_in);
        particles.push_back(p_in);
    }

    // External electric field at point r=(x,y,z)
    vec external_E_field(int i)
    {
        vec P = vec(3); //vector for positions

        P(0) = -1.; // positions for Electric-field
        P(1) = -1.;
        P(2) = 2.;

        vec r = particles[i].r;
        vec electricfield = (((-1 * V0)) * P) % r; //TODO d^2 not sure

        return electricfield;
    }

    // External magnetic field at point r=(x,y,z)
    vec external_B_field(int i)
    {
        vec magneticfield = vec(3).fill(0.);
        magneticfield(2) = B0;
        return magneticfield;
    }

    // coulombic interactions Force on particle_i from particle_j
    vec force_particle(int i, int j)
    {
        Particle pi = particles[i];
        Particle pj = particles[j];

        vec r = pi.r - pj.r;

        vec r3 = abs(r) % abs(r) % abs(r);

        vec force = ke * (pi.q * pj.q) / r3 % r;

        return force;
    }

    vec coulomb_interaction(int i, int j)
    {
        if (i != j)
        {
            return force_particle(i, j);
        }
        else
        {
            return vec(3).fill(0);
        }
    }

    //total force onparticle i due to external fields
    vec total_force_external(int i)
    {

        vec F = vec(3).fill(0);
        vec E = external_E_field(i);
        vec B = external_B_field(i);

        Particle p = particles[i];
        vec v = p.v;
        double q = p.q;

        F = q * E + cross(q * v, B);
        return F;
    }
    
    //force on particle i from other particles
    vec total_force_particles(int i)
    {
        vec F = vec(3).fill(0);
        for (int j = 0; j < n; j++)
        {
             F += coulomb_interaction(i, j);
        }
        return F;
    }

    //force on particle due to fields and particles
    vec total_force(int i)
    {
        vec F = total_force_particles(i) + total_force_external(i);

        return F;
    }

    void evolve_RK4(double dt, mat& R, mat& V)
    {

        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];

            vec F = total_force(i);
            vec acceleration = F / p.m;

            vec k0_velocity = p.v; // particle initial velocity
            vec k0_position = p.r; // particle initial position

            // K1 velocity and position
            vec k1_v = acceleration * dt;
            vec k1_velocity = k0_velocity + k1_v / 2;

            vec k1_r = k1_velocity * dt;
            vec k1_position = k0_position + k1_r / 2;

            // K2
            vec k2_v = k0_velocity + k1_v / 2;
            vec k2_velocity = k1_velocity + k2_v / 2;

            vec k2_r = k0_position + k1_r / 2;
            vec k2_position = k1_position + k2_r / 2;

            // K3
            vec k3_v = k1_velocity + k2_v / 2;
            vec k3_velocity = k2_velocity + k3_v / 2;

            vec k3_r = k1_position + k2_r / 2;
            vec k3_position = k1_position + k3_r / 2;

            // K4
            vec k4_v = k2_velocity + k3_v;
            vec k4_velocity = k3_velocity + k4_v;

            vec k4_r = k2_position + k3_r;
            vec k4_position = k3_position + k4_r;

            V.col(i) = k4_velocity + (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6;
            R.col(i) = k4_position + (k1_r + 2 * k2_r + 2 * k3_r + k4_r) / 6;
        }

    }

    void evolve_forward_Euler(double dt, mat& R, mat& V)
    {
        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];
            vec F = total_force(i);

            vec a = F / p.m;

            V.col(i) = p.v + a * dt;
            R.col(i) = p.r + p.v * dt;
        }
    }
};
