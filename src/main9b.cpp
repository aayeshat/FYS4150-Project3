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
    double t = 100.;
    double N = 100000.;
    double dt = t / N;

    double B0 = 96.5;
    double V0 = 9.65; //v0/d^2
    double d = 10e4;
    int number_of_particles = 2;

    PenningTrap trap(B0, V0, d, number_of_particles);
    trap.interaction = true;

    for (int i = 0; i < number_of_particles; i++)
    {
        vec r = vec(3, fill::randu) * 1.0e4;
        vec v = vec(3, fill::randu);
        Particle particle_i(1., 40.078, r, v);
        trap.add_particle(particle_i);
    }

    mat position1_evol_inter(N, 3);
    mat position2_evol_inter(N, 3);
    mat position1_evol_nointer(N, 3);
    mat position2_evol_nointer(N, 3);

    mat r_step(3, number_of_particles);
    mat v_step(3, number_of_particles);

    // with coloumbib interactions
    cube R(3, N, trap.particles.size(), fill::zeros);

    ofstream out;
    out.open("./out/xy_inter_1_2.txt");
    for (int k = 0; k < N; k++)
    {
        out << setprecision(4) << scientific << (k * dt);
        for (int i = 0; i < trap.particles.size(); i++)
        {
            if (k == 0)
            {
                R.slice(i).col(k) = trap.particles[i].r;

                r_step.col(i) = trap.particles[i].r;
                v_step.col(i) = trap.particles[i].v;
            }
            else
            {
                trap.evolve_RK4(dt, i, r_step, v_step);
                R.slice(i).col(k) = r_step.col(i);
            }

            mat col = R.slice(i).col(k).t();
            for (int j = 0; j < 3; j++)
            {
                out << "   " << col(j);
            }
        }

        out << endl;
    }

    out.close();

    return 0;
}