#include "PenningTrap.hpp"

#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, int n_in)
{
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
    n = n_in;
}

void PenningTrap::add_particle(Particle p_in)
{
    particles.push_back(p_in);
}

vec PenningTrap::external_E_field(int i, double t)
{
    Particle p = particles[i];
    vec r = p.r;
    vec E = vec(3).fill(0);

    if (norm(r) <= d)
    {
        E[0] = V0 * (1 + f * cos(omega_v * t)) * r[0] / pow(d, 2);
        E[1] = V0 * (1 + f * cos(omega_v * t)) * r[1] / pow(d, 2);
        E[2] = -V0 * (1 + f * cos(omega_v * t)) * 2 * r[2] / pow(d, 2);
    }
    return E;
}

// External magnetic field at point r=(x,y,z)
vec PenningTrap::external_B_field(int i)
{
    vec magneticfield = vec(3).fill(0.);
    magneticfield(2) = B0;
    return magneticfield;
}

// coulombic interactions Force on particle_i from particle_j
vec PenningTrap::force_particle(int i, int j)
{
    Particle pi = particles[i];
    Particle pj = particles[j];

    vec r = pi.r - pj.r;

    vec r3 = abs(r) % abs(r) % abs(r);

    vec force = ke * (pi.q * pj.q) / r3;

    return force;
}

//total force on particle i due to external fields

vec PenningTrap::total_force_external(int i, double t)
{
    Particle p = particles[i];
    vec v = p.v;
    double q = p.q;

    vec F = vec(3).fill(0);

    vec E = external_E_field(i, t);

    vec B = external_B_field(i);

    F(0) = q * (E(0) + v(1) * B(2));
    F(1) = q * (E(1) - v(0) * B(2));
    F(2) = q * E(2);

    return F;
}

//force on particle i from other particles
vec PenningTrap::total_force_particles(int i)
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
vec PenningTrap::total_force(int i, double t)
{
    vec F = vec(3).fill(0);
    if (interaction)
    {
        for (int k = 0; k <= 2; k++)
        {
            F(k) = total_force_particles(i)(k) + total_force_external(i, t)(k);
        }
    }
    else
    {
        for (int k = 0; k <= 2; k++)
        {
            F(k) = total_force_external(i, t)(k);
        }
    }
    return F;
}

void PenningTrap::evolve_RK4(double dt, double t)
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
        k1_v.col(i) = dt * total_force(i, t + 0.5 * dt) / p.m;
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
        k2_v.col(i) = dt * total_force(i, t + 0.5 * dt) / p.m;

        particles[i].r = initial_r.col(i) + 0.5 * k2_r.col(i);
        particles[i].v = initial_v.col(i) + 0.5 * k2_v.col(i);
    }

    //K3
    for (int i = 0; i < n; i++)
    {
        Particle p = particles[i];

        k3_r.col(i) = dt * p.v;
        k3_v.col(i) = dt * total_force(i, t + 0.5 * dt) / p.m;
        particles[i].r = initial_r.col(i) + k3_r.col(i) * 0.5;
        particles[i].v = initial_v.col(i) + k3_v.col(i) * 0.5;
    }

    //K4
    for (int i = 0; i < n; i++)
    {
        Particle p = particles[i];
        k4_r.col(i) = dt * p.v;
        k4_v.col(i) = dt * total_force(i, t + 0.5 * dt) / p.m;

        particles[i].r = initial_r.col(i) + (k1_r.col(i) + 2 * k2_r.col(i) + 2 * k3_r.col(i) + k4_r.col(i)) / 6;
        particles[i].v = initial_v.col(i) + (k1_v.col(i) + 2 * k2_v.col(i) + 2 * k3_v.col(i) + k4_v.col(i)) / 6;
    }
}


void PenningTrap::evolve_forward_Euler(double dt, double t)
{
    for (int i = 0; i < n; i++)
    {
        Particle p = particles[i];
        vec F = total_force(i, t);
        vec a = F / p.m;

        particles[i].v = p.v + a * dt;
        particles[i].r = p.r + p.v * dt;
    }
}
