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
    double tf = 100.;
    double N = 1000.;
    double dt = tf / N;
    vec t = linspace(0,tf,dt);
    double B0 = 96.5;
    double V0 = 9.65; //v0/d^2
    double d = 10e4;

    int number_of_particles = 10;

    PenningTrap trap(B0, V0, d, number_of_particles);
    trap.interaction = true; //switch for interaction true (for interactions) or false (without coulombic interactions)

    vec r = vec(3).randn() * 0.1 * trap.d; // random initial position
    vec v = vec(3).randn() * 0.1 * trap.d; // random initial velocity

    Particle particle_i(1., 40.078, r, v);
    trap.add_particle(particle_i);
}


/* 
if (trap.interaction)
{
    position_out_filename = "./out/r_xy_inter_1_2_initial_condition.txt";
    velocity_out_filename = "./out/v_xy_inter_1_2_initial_condition.txt";
}
else
{
    position_out_filename = "./out/r_xy_nointer_1_2_initial_condition.txt";
    velocity_out_filename = "./out/v_xy_nointer_1_2_initial_condition.txt";
}
position_out.open(position_out_filename);
velocity_out.open(velocity_out_filename);
 */

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

//position_out.close();
//velocity_out.close();

return 0;
}