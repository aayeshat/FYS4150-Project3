#include "Particle.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

class PenningTrap
{
public:
    double B0;
    double V0;
    double d;
    int n;
    double ke = 1.38935333e5;
    vector<Particle> particles;

    PenningTrap(double B0_in, double V0_in, double d_in, int n_in) // Constructor
    {
        B0 = B0_in; //magnetic filed
        V0 = V0_in; //Applied Potential
        d = d_in;   //dimension
        n = n_in;
    }
    void add_particle(Particle p_in)
    {
        //particles.insert(particles.begin(), p_in);
        particles.push_back(p_in);
    }

    // External electric field at point r=(x,y,z)
    vec external_E_field(int i)
    {
        vec P = vec(3); //vector for positions

        P(0) = -1.; // positions for Electric-field
        P(1) = -1.;
        P(2) = 2.;

        vec r = particles[i].r;
        vec electricfield = (((-1 * V0) / pow(d, 2)) * P) % r; //TODO d^2 not sure

        return electricfield;
    }

    // External magnetic field at point r=(x,y,z)
    vec external_B_field(int i)
    {
        vec magneticfield = vec(3).fill(0.);
        magneticfield(2) = B0;
        return magneticfield;
    }

    // Force on particle_i from particle_j
    vec force_particle(int i, int j)
    {
        Particle pi = particles[i];
        Particle pj = particles[j];

        vec r = pi.r - pj.r;

        vec r3 = abs(r) % abs(r) % abs(r);

        vec force = ke * (pi.q * pj.q) / r3 % r;

        return force;
    }

    vec total_force_external(int i)
    {

        vec F = vec(3).fill(0);
        vec E = external_E_field(i);
        vec B = external_B_field(i);

        Particle p = particles[i];
        vec v = p.v;
        double q = p.q;

        F = q * E + cross(q * v, B);
        return F;
    }

    vec total_force_particles(int i)
    {
        vec F = vec(3).fill(0);
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                F += force_particle(i, j);
            }
        }
        return F;
    }

    vec total_force(int i)
    {
        vec F = total_force_particles(i) + total_force_external(i);

        return F;
    }

    void evolve_RK4(double dt)
    {
    }

    void evolve_forward_Euler(double dt)
    {
        mat R = mat(3, n).fill(0);
        mat V = mat(3, n).fill(0);

        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];
            vec F = total_force(i);

            vec a = F / p.m;

            V.col(i) = p.v + a * dt;
            R.col(i) = p.r + p.v * dt;
        }

        V.print("V= ");
        R.print("R= ");

        int width = 12;
        int prec = 4;
        string filename = "./out/evolve_forward_Euler.txt";
        ofstream ofile;
        ofile.open(filename);

        cout << setw(width) << setprecision(prec) << scientific << "R"
             << setw(width) << setprecision(prec) << scientific << "V"
             << endl;

        for (int i = 0; i < V.size(); i++)
        {

            cout << setw(width) << setprecision(prec) << scientific << R(i)
                 << setw(width) << setprecision(prec) << scientific << V(i)
                 << endl;

            ofile << setw(width) << setprecision(prec) << scientific << R(i)
                  << setw(width) << setprecision(prec) << scientific << V(i)
                  << endl;
        }
        ofile.close();
    }
};
