#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[])
{
    double q = 20;                // charge of Ca+ particle [e]
    double m = 40.078;            // atomic mass of Ca+ [u]

    double B0 = 9.65e1;           // magnetic field strength
    double V0 = 9.65e8;           // applied potential
    double d = 10e4;              // characteristic dimension

    double ke = 1.38935333e5;     // Couloumb constant

    int n = 2;                    // number of particles
    int dim = 3;                  // dimension (x,y,z)

    
    PenningTrap penningtrap = PenningTrap(B0, V0, d);    // calling penningtrap

   
    return 0;

}