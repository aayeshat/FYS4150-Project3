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

    mat r_step;
    mat v_step;
    bool interaction = true;

    PenningTrap(double B0_in, double V0_in, double d_in, int n_in) // Constructor
    {
        B0 = B0_in; //magnetic filed
        V0 = V0_in; //Applied Potential
        d = d_in;   //dimension
        n = n_in;

        r_step = mat(3, n);
        v_step = mat(3, n);
    }

    void add_particle(Particle p_in)
    {
        //particles.insert(particles.begin(), p_in);
        particles.push_back(p_in);
    }

    // External electric field at point r=(x,y,z)
    vec external_E_field(int i)
    {
        Particle p = particles[i];
        vec r = p.r;

        double V_d = 9.65;
        vec E = vec(3).fill(0);

        E(0) = V_d * r(0);
        E(1) = V_d * r(1);
        E(2) = -(2 * V_d) * r(2);

        return E;
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
        int ke = 1.38935333e5;
        Particle pi = particles[i];
        Particle pj = particles[j];

        vec r = pi.r - pj.r;

        vec r3 = abs(r) % abs(r) % abs(r);

        vec force = ke * (pi.q * pj.q) / r3;

        return force;
    }

    //total force on particle i due to external fields

    vec total_force_external(int i)
    {
        Particle p = particles[i];
        vec v = p.v;
        double q = p.q;

        vec F = vec(3).fill(0);

        vec E = external_E_field(i);

        vec B = external_B_field(i);
        // F = q * E + cross(q * v, B);

        F(0) = q * (E(0) + v(1) * B(2));
        F(1) = q * (E(1) - v(0) * B(2));
        F(2) = q * E(2);

        return F;
    }

    //force on particle i from other particles
    vec total_force_particles(int i)
    {
        vec F = vec(3).fill(0);
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                return force_particle(i, j);
            }
        }

        return F;
    }

    //force on particle due to fields and particles
    vec total_force(int i)
    {
        vec F = vec(3).fill(0);
        if (interaction)
        {
            for (int k = 0; k <= 2; k++)
            {
                F(k) = total_force_particles(i)(k) + total_force_external(i)(k);
            }
        }
        else
        {
            for (int k = 0; k <= 2; k++)
            {
                F(k) = total_force_external(i)(k);
            }
        }
        return F;
    }

    void evolve_RK4(double dt)
    {
        mat initial_r(3, n);
        mat initial_v(3, n);

        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];
            initial_r.col(i) = p.r; // particle initial position
            initial_v.col(i) = p.v; // particle initial velocity
        }

        mat k1_r(3, n);
        mat k2_r(3, n);
        mat k3_r(3, n);
        mat k4_r(3, n);

        mat k1_v(3, n);
        mat k2_v(3, n);
        mat k3_v(3, n);
        mat k4_v(3, n);

        //K1
        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];

            k1_r.col(i) = dt * p.v;
            k1_v.col(i) = dt * total_force(i) / p.m;
        }

        for (int i = 0; i < n; i++)
        {
            particles[i].r = initial_r.col(i) + k1_r.col(i) / 2;
            particles[i].v = initial_v.col(i) + k1_v.col(i) / 2;
        }

        //K2
        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];

            k2_r.col(i) = dt * p.v;
            k2_v.col(i) = dt * total_force(i) / p.m;

            particles[i].r = initial_r.col(i) + 0.5 * k2_r.col(i);
            particles[i].v = initial_v.col(i) + 0.5 * k2_v.col(i);
        }

        //K3
        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];

            k3_r.col(i) = dt * p.v;
            k3_v.col(i) = dt * total_force(i) / p.m;
            particles[i].r = initial_r.col(i) + k3_r.col(i) * 0.5;
            particles[i].v = initial_v.col(i) + k3_v.col(i) * 0.5;
        }

        //K4
        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];
            k4_r.col(i) = dt * p.v;
            k4_v.col(i) = dt * total_force(i) / p.m;

            particles[i].r = initial_r.col(i) + (k1_r.col(i) + 2 * k2_r.col(i) + 2 * k3_r.col(i) + k4_r.col(i)) / 6;
            particles[i].v = initial_v.col(i) + (k1_v.col(i) + 2 * k2_v.col(i) + 2 * k3_v.col(i) + k4_v.col(i)) / 6;
        }

        for (int i = 0; i < n; i++)
        {
            r_step.col(i) = particles[i].r;
            v_step.col(i) = particles[i].v;
        }
    }

    void evolve_forward_Euler(double dt)
    {

        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];
            vec F = total_force(i);
            vec a = F / p.m;

            particles[i].v = p.v + a * dt;
            particles[i].r = p.r + p.v * dt;
        }

        for (int i = 0; i < n; i++)
        {
            r_step.col(i) = particles[i].r;
            v_step.col(i) = particles[i].v;
        }
    }
};
