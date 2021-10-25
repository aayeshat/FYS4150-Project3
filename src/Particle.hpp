#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

class Particle
{
public:
    double q; //charge
    double m; //mass
    vec r; //position
    vec v; //velocity

    Particle(double q_in, double m_in, vec r_in, vec v_in) //constructor
    {
        q = q_in;
        m = m_in;
        r = r_in;
        v = v_in;
    }

    void print_particle() // funtion
    {
        cout << "charge = " << q << endl;
        cout << "mass = " << m << endl;
        r.print("position = ");
        v.print("velocity = ");
    }
};