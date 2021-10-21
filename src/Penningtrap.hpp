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

    mat R;
    mat V;
    bool interaction = true;

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
        Particle pi = particles[i];
        Particle pj = particles[j];

        vec r = pi.r - pj.r;

        vec r3 = abs(r) % abs(r) % abs(r);

        vec force = ke * (pi.q * pj.q) / r3 % r;

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

    void evolve_RK4(double dt, int i, mat &R, mat &V)
    {
        vec position_;
        vec velocity_;
        int number_of_particles;
        int k;

        Particle p_old = particles[i];

        double m = p_old.m;
/*         mat initial_r = mat(3,number_of_particles);
        mat initial_v = mat(3,number_of_particles);

        initial_r.col(i) = p_old.r;
        initial_v.col(i) = p_old.v; */

        mat K_r = mat(3, number_of_particles);
        mat K_v = mat(3, number_of_particles);
        mat K_sum_r = mat(3, number_of_particles);
        mat K_sum_v = mat(3, number_of_particles);

        for (int k = 0; k <= number_of_particles; k++)
        { //K1

            Particle p_i = particles[i];
            Particle p_old_i = p_old;

            K_r.col(i) = dt * p_i.v;
            K_v.col(i) = dt * total_force(i) / m;

            p_i.r = p_old_i.r + K_r.col(i) / 2;
            p_i.v = p_old_i.v + K_v.col(i) / 2;
        }

        K_sum_r += K_r;
        K_sum_v += K_v;

        for (int k = 0; k <= number_of_particles; k++)
        { //K2

            Particle p_i = particles[i];
            Particle p_old_i = p_old;

            K_r.col(i) = dt * p_i.v;
            K_v.col(i) = dt * total_force(i) / m;

            p_i.r = p_old_i.r + K_r.col(i) / 2;
            p_i.v = p_old_i.v + K_v.col(i) / 2;
        }


        K_sum_r += 2 * K_r;
        K_sum_v += 2 * K_v;

        //K3
        for (int k = 0; k <= number_of_particles; k++)
        {

            Particle p_i = particles[i];
            Particle p_old_i = p_old;

            K_r.col(i) = dt * p_i.v;
            K_v.col(i) = dt * total_force(i) / m;

            p_i.r = p_old_i.r + K_r.col(i) / 2;
            p_i.v = p_old_i.v + K_v.col(i) / 2;
        }

        p_old.r = particles[i].r;
        p_old.v = particles[i].v;

        K_sum_r += 2 * K_r;
        K_sum_v += 2 * K_v;

        //K4
        for (int k = 0; k <= number_of_particles; k++)
        {
            Particle p_i = particles[i];
            Particle p_old_i = p_old;

            K_r.col(i) = dt * p_i.v;
            K_v.col(i) = dt * total_force(i) / m;

            p_i.r = p_old_i.r + K_r.col(i);
            p_i.v = p_old_i.v + K_v.col(i);
        }

        K_sum_r += K_r;
        K_sum_v += K_v;

        vec r_step = p_old.r + K_sum_r / 6;
        vec v_step = p_old.v + K_sum_v / 6;

        R.col(i) = r_step;
        V.col(i) = v_step;
    }

    void evolve_RK4(double dt)
    {
        R = mat(3, n).fill(0);
        V = mat(3, n).fill(0);

        for (int i = 0; i < n; i++)
        {

            Particle p = particles[i];
            vec initial_r = p.r; // particle initial position
            vec initial_v = p.v; // particle initial velocity

            evolve_RK4(dt, i, R, V);

            // reset to initial
            particles[i].r = initial_r;
            particles[i].v = initial_v;
        }
    }

    /* void evolve_RK4(double dt, int i, mat &R, mat &V)
{
vec position_;
vec velocity_;

Particle p = particles[i];
double m = p.m;

vec initial_r = p.r; // particle initial position
vec initial_v = p.v; // particle initial velocity

//K1
vec k1_r = dt * p.v;
vec k1_v = dt * total_force(i) / m;
particles[i].r = initial_r + k1_r / 2;
particles[i].v = initial_v + k1_v / 2;

//K2
vec k2_r = dt * particles[i].v;
vec k2_v = dt * total_force(i) / m;

particles[i].r = initial_r + k2_r / 2;
particles[i].v = initial_v + k2_v / 2;

//K3
vec k3_r = dt * particles[i].v;
vec k3_v = dt * total_force(i) / m;

particles[i].r = initial_r + k3_r / 2;
particles[i].v = initial_v + k3_v / 2;

//K4
vec k4_r = dt * particles[i].v;
vec k4_v = dt * total_force(i) / m;

vec r_step = initial_r + (k1_r + 2 * k2_r + 2 * k3_r + k4_r) / 6;
vec v_step = initial_v + (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6;

R.col(i) = r_step;
V.col(i) = v_step;
}

void evolve_RK4(double dt)
{
R = mat(3, n).fill(0);
V = mat(3, n).fill(0);

for (int i = 0; i < n; i++)
{

    Particle p = particles[i];
    vec initial_r = p.r; // particle initial position
    vec initial_v = p.v; // particle initial velocity

    evolve_RK4(dt, i, R, V);

    // reset to initial
    particles[i].r = initial_r;
    particles[i].v = initial_v;
}
}
 */
    void evolve_forward_Euler(double dt, int i, mat &R, mat &V)
    {
        Particle p = particles[i];
        vec F = total_force(i);

        vec a = F / p.m;

        V.col(i) = p.v + a * dt;
        R.col(i) = p.r + p.v * dt;
    }

    void evolve_forward_Euler(double dt)
    {
        R = mat(3, n).fill(0);
        V = mat(3, n).fill(0);
        for (int i = 0; i < n; i++)
        {
            evolve_forward_Euler(dt, i, R, V);
        }
    }
};
