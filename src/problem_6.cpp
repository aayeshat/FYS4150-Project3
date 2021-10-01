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

    Particle(double charge_in, double mass_in, vec position_in, vec velocity_in) //constructer
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

int main()
{

    Particle p(-1., -2., vec(2, fill::eye), vec(3, fill::randu));
    
    p.print_particle();

    cout << "charge inside class = " << p.charge << endl;

    return 0;
}