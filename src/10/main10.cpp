#include "Penningtrap.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    int t0 = 0;
    double t = 500.;
    double N = 500.;
    double dt = t / N;
    double B0 = 9.65e1;
    double V0 = 0.0025 * 9.64852558 * 1.0e7;
    double d = 500;
    int number_of_particles = 100;

    PenningTrap trap(B0, V0, d, number_of_particles);

    trap.interaction = false; //switch for interaction true (for interactions) or false (without coulombic interactions)

    arma_rng::set_seed_random();

    for (int i = 0; i < number_of_particles; i++)
    {
        vec r = vec(3).randn() * 0.1 * trap.d;
        vec v = vec(3).randn() * 0.1 * trap.d;

        Particle particle_i(1., 40.078, r, v);
        trap.add_particle(particle_i);
    }

    string out_filename;

    if (trap.interaction)
    {
        out_filename = "./out/10_inter__.txt";
    }
    else
    {
        out_filename = "./out/10_nointer__.txt";
    }

    ofstream out;
    out.open(out_filename);

    for (double omega_v = 0.2; omega_v <= 2.5; omega_v += 0.02)
    {
        trap.omega_v = omega_v;
        for (int k = 0; k < N; k++)
        {
            trap.evolve_RK4(dt, k + dt);
        }

        int number_inside = 0;
        for (int i = 0; i < number_of_particles; i++)
        {
            vec r = trap.particles[i].r;
            cout << norm(r) << " < " << trap.d << endl;
            cout << "-------------" << endl;
            if (norm(r) < trap.d)
            {
                number_inside += 1;
            }
        }

        out << scientific << omega_v << "   " << fixed << number_inside << endl;

        for (int i = 0; i < number_of_particles; i++)
        {
            trap.particles[i].r = vec(3).randn() * 0.1 * trap.d;
            trap.particles[i].v = vec(3).randn() * 0.1 * trap.d;
        }
    }

    out.close();

    return 0;
}