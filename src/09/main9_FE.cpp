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
    double N = 1000.; //stepsize is 0.1
    double dt = t / N;

    double B0 = 96.5;
    double V0 = 9.65e8; 
    double d = 10e4;
    int number_of_particles = 1;

    PenningTrap trap(B0, V0, d, number_of_particles);
    
    trap.interaction = false; //switch for interaction true (for interactions) or false (without coulombic interactions)
    

    for (int i = 0; i < number_of_particles; i++)
    {
        vec r, v;
        if (i == 0)
        {
            r = vec(3).fill(0)*10e4;

            r(2) = 10.;
            r(0) = 10.;

            v = vec(3).fill(0);
            v(1) = 10.;
        }
        else
        {
            r = vec(3, fill::randn);
            v = vec(3).randn()*0.001;
        }

        Particle particle_i(1., 40.078, r, v);
        trap.add_particle(particle_i);
    }

    ofstream position_out;
    string position_out_filename;

    if (trap.interaction)
    {
        position_out_filename = "./out/data/r_1_euler(0.1).txt";
    }
    else
    {
        position_out_filename = "./out/data/r_1_euler(0.1).txt";
    }
    position_out.open(position_out_filename);

    for (int k = 0; k < N; k++)
    {

        mat r_step;
        mat v_step;

        position_out << setprecision(4) << scientific << (k * dt);

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
            trap.evolve_forward_Euler(dt);
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
            }
        }

        position_out << endl;
    }

    position_out.close();

    return 0;
}