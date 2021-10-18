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
    double N = 10000.;
    double dt = t / N;

    double B0 = 96.5;
    double V0 = 9.65; //v0/d^2
    double d = 10e4;
    int number_of_particles = 1;

    PenningTrap trap(B0, V0, d, number_of_particles);

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

    // with interactions
    cube R(3, N, trap.particles.size(), fill::zeros);

    ofstream myfile_t;
    myfile_t.open("./out/time.txt");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < trap.particles.size(); j++)
        {
            if (i == 0)
            {
                R.slice(j).col(i) = trap.particles[j].r;
                r_step.col(j) = trap.particles[j].r;
                v_step.col(j) = trap.particles[j].v;
            }
            else
            {
                trap.evolve_RK4(dt, j, r_step, v_step);
                R.slice(j).col(i) = r_step.col(j);
            }
        }

        myfile_t << (i * dt) << endl;
    }


    ofstream myfile1_inter;
    myfile1_inter.open("./out/r_values.txt");
    myfile1_inter << R.slice(0).t();
    myfile1_inter.close();


    return 0;
}