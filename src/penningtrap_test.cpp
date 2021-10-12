#include "Penningtrap.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{

    double B0 = 9.65e1;
    double V0 = 9.65e8;
    double d = 10e4;
    int n = 2;

    PenningTrap pt(B0, V0, d, n);

    Particle particle1(1., 40.078, vec(3, fill::randu), vec(3, fill::randu)); //TODO random position and velocity for now
    pt.add_particle(particle1);

    Particle particle2(1., 40.078, vec(3, fill::randu), vec(3, fill::randu)); //TODO  random position and velocity for now
    pt.add_particle(particle2);

    vec electricfield = pt.external_E_field(0);
    electricfield.print("p1 electricfield");
    vec magneticfield = pt.external_B_field(0);
    magneticfield.print("p1 magneticfield");

    vec force = pt.force_particle(0, 1);
    force.print("force");

    cout << endl
         << endl
         << "evolve_RK4 "
         << endl
         << endl;
    pt.evolve_RK4(.1);

    cout << endl
         << endl
         << "evolve_forward_Euler "
         << endl
         << endl;
    pt.evolve_forward_Euler(.1);

    return 0;
}