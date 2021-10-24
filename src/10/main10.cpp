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
    double N = 1000.;
    double dt = t / N;

    double B0 = 96.5;
    double V0 = 0.0025*9.64852558; //v0/d^2
    double d = 500;
    int number_of_particles = 10;

    PenningTrap trap(B0, V0, d, number_of_particles);

    trap.interaction = true; //switch for interaction true (for interactions) or false (without coulombic interactions)

    for (int i = 0; i < number_of_particles; i++)
    {
         vec  r =  vec(3).randn() * 0.1 * trap.d;
         vec  v =  vec(3).randn() * 0.1 * trap.d;

        Particle particle_i(1., 40.078, r, v);
        trap.add_particle(particle_i);
    }

    string out_filename;

    if (trap.interaction)
    {
        out_filename = "./out/10_inter__0.1.txt";
    }
    else
    {
        out_filename = "./out/10_nointer__0.txt";
    }

    ofstream out;
    out.open(out_filename);

    for (double omega_v = 0.2; omega_v <= 2.5; omega_v += 0.02)
    {
        trap.omega_v = omega_v;
        for (int k = 0; k < N; k++)
        {
            trap.evolve_RK4(dt, dt * k);
        }

        double counter = 0;
        for (int i = 0; i < number_of_particles; i++)
        {
            mat r_step = trap.r_step / number_of_particles;
            vec r = trap.r_step.col(i);
            double sq = sqrt(pow(r(0), 2) + pow(r(1), 2) + pow(r(2), 2));
            if (sq < trap.d)
            {
                counter += 1;
            }
        }
        out << setprecision(4) << scientific << omega_v << "   " << counter << endl;
    }

    out.close();

    return 0;
}