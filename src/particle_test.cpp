#include "Particle.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    Particle p(1., 40., vec(3, fill::ones), vec(3, fill::randu));

    p.print_particle();

    cout << "charge inside class = " << p.q << endl;

    return 0;
}