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
    int t = 100;
    int N = 1002;

    double B0 = 96.5;
    double V0 = 9.65; //v0/d^2
    double d = 10e4;
    int number_of_particles = 2;

    PenningTrap trap(B0, V0, d, number_of_particles);

    for (int i = 0; i < number_of_particles; i++)
    {
        Particle particle_i(1., 40.078, vec(3, fill::randu), vec(3, fill::randu));
        trap.add_particle(particle_i);
    }

    vec time = linspace(0, t, N);

    int width = 12;
    int prec = 4;

    ofstream ofileV;
    ofileV.open("./out/data_V2.txt");

    ofstream ofileR;
    ofileR.open("./out/data_R2.txt");

    for (int k = 0; k < N; k++)
    {
        double dt = time(k);
     
        trap.evolve_RK4(dt);
        //trap.evolve_forward_Euler(dt);

        ofileV << setw(width) << setprecision(prec) << scientific << dt;
        for (int i = 0; i < trap.V.size(); i++)
        {
            ofileV << setw(width) << setprecision(prec) << scientific << trap.V(i);
        }
        ofileV << endl;

        ofileR << setw(width) << setprecision(prec) << scientific << dt;
        for (int i = 0; i < trap.R.size(); i++)
        {
            ofileR << setw(width) << setprecision(prec) << scientific <<trap. R(i);
        }
        ofileR << endl;
    }

    ofileV.close();
    ofileR.close();
    return 0;
}