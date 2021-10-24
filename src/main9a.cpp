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
    double t = 1000.;
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
    
    ofstream position_out;
    position_out.open("./out/r_values.txt");
    ofstream velocity_out;
    velocity_out.open("./out/v_values.txt");

    for (int k = 0; k < N; k++)
    {
        mat r_step;
        mat v_step;

        position_out << setprecision(4) << scientific << (k * dt);
        velocity_out << setprecision(4) << scientific << (k * dt);

        if (k == 0)
        {
            r_step = mat(3, number_of_particles);
            v_step = mat(3, number_of_particles);

            for (int i = 0; i < number_of_particles; i++)
            {
                r_step.col(i) = trap.particles[i].r;
                v_step.col(i) = trap.particles[i].v;
            }
        }
        else
        {
            trap.evolve_RK4(dt);
            r_step = trap.r_step;
            v_step = trap.v_step;
        }

        for (int i = 0; i < number_of_particles; i++)
        {
            vec p_r = r_step.col(i);
            vec p_v = v_step.col(i);
            for (int j = 0; j < 3; j++)
            {
                position_out << "   " << p_r(j);
                velocity_out << "   " << p_v(j);
            }
        }

        position_out << endl;
        velocity_out << endl;
    }

    position_out.close();
    velocity_out.close();
    return 0;
}