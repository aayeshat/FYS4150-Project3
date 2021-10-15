#include "Penningtrap.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{

    int t = 100;
    int N = 1000;
    double dt = t * (1. / N);

    double B0 = 9.65e1;
    double V0 = 9.65e8;
    double d = 10e4;
    int number_of_particles = 1;

    PenningTrap trap(B0, V0, d, number_of_particles);

    for (int i = 0; i < number_of_particles; i++)
    {
        Particle particle_i(1., 40.078, vec(3, fill::randu), vec(3, fill::randu));
        trap.add_particle(particle_i);
    }

    vec time = linspace(0, t, t / dt);

    ofstream ofile;
    ofile.open("./out/data.txt");
    int width = 12;
    int prec = 4;

    cout << setw(width) << setprecision(prec) << scientific << "dt"
         << setw(width) << setprecision(prec) << scientific << "R"
         << setw(width) << setprecision(prec) << scientific << "V"
         << endl;

    ofile <<"# "<< setprecision(prec) << scientific << "dt"
         << setw(width) << setprecision(prec) << scientific << "R"
         << setw(width) << setprecision(prec) << scientific << "V"
         << endl;

    for (int k = 0; k < N; k++)
    {
        double dt1 = time(k);
        mat R = mat(3, number_of_particles).fill(0);
        mat V = mat(3, number_of_particles).fill(0);
        trap.evolve_RK4(dt1, R, V);

        for (int i = 0; i < R.size(); i++)
        {

            cout << setw(width) << setprecision(prec) << scientific << dt1
                 << setw(width) << setprecision(prec) << scientific << R(i)
                 << setw(width) << setprecision(prec) << scientific << V(i)
                 << endl;

            ofile << setw(width) << setprecision(prec) << scientific << dt1
                  << setw(width) << setprecision(prec) << scientific << R(i)
                  << setw(width) << setprecision(prec) << scientific << V(i)
                  << endl;
        }

        cout << "---------------------------" << endl;
    }

    ofile.close();
    return 0;
}