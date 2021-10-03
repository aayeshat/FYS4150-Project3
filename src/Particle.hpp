#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

class Particle
{
public:
    double charge;
    double mass;
    vec position;
    vec velocity;

    Particle(double charge_in, double mass_in, vec position_in, vec velocity_in) //constructor
    {
        charge = charge_in;
        mass = mass_in;
        position = position_in;
        velocity = velocity_in;
    }

    void print_particle() // funtion
    {
        cout << "charge = " << charge << endl;
        cout << "mass = " << mass << endl;
        position.print("position = ");
        velocity.print("velocity = ");
    }
};