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

    Particle p1(1., 40.078, vec(3, fill::randu), vec(3, fill::randu)); //TODO random position and velocity for now
    pt.add_particle(p1);

    Particle p2(1., 40.078, vec(3, fill::randu), vec(3, fill::randu)); //TODO  random position and velocity for now
    pt.add_particle(p2);

    vec electricfield = pt.external_E_field(0);
    electricfield.print("p1 electricfield");
    vec magneticfield = pt.external_B_field(0);
    magneticfield.print("p1 magneticfield");

    vec electricfield1 = pt.external_E_field(1);
    electricfield1.print("p2 electricfield");
    vec magneticfield1 = pt.external_B_field(1);
    magneticfield1.print("p2 magneticfield");

    vec force = pt.force_particle(0, 1);
    force.print("force");

    return 0;
}