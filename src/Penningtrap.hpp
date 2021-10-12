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
        mat R = mat(3, n).fill(0);
        mat V = mat(3, n).fill(0);
        for (int i = 0; i < n; i++)
        {
            Particle p = particles[i];

            vec F = total_force(i);
            vec acceleration = F / p.m;

            vec velocity = p.v; // particle initial velocity
            vec k0_velocity = velocity;

            vec position = p.r; // particle initial position
            vec k0_position = position;

            // K1 velocity and position
            vec k1_v = acceleration * dt;
            velocity = velocity + k1_v / 2;
            vec k1_velocity = velocity;

            vec k1_r = velocity * dt;
            position = position + k1_r / 2;
            vec k1_position = position;

            // K2
            vec k2_v = k0_velocity + k1_v / 2;
            velocity = velocity + k2_v / 2;
            vec k2_velocity = velocity;

            vec k2_r = k0_position + k1_r / 2;
            position = position + k2_r / 2;
            vec k2_position = position;

            // K3 
            vec k3_v = k1_velocity + k2_v / 2;
            velocity = velocity + k3_v / 2;
            vec k3_r = k1_position + k2_r / 2;
            position = position + k3_r / 2;

            // K4 
            vec k4_v = k2_velocity + k3_v / 2;
            velocity = velocity + k4_v / 2;
            vec k4_r = k2_position + k3_r / 2;
            position = position + k4_r / 2;


            V.col(i) = velocity + (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6;
            R.col(i) = position + (k1_r + 2 * k2_r + 2 * k3_r + k4_r) / 6;
        }

        V.print("evolve_RK4 V= ");
        R.print("evolve_RK4 R= ");

        printAndSaveToFile(V, R, "./out/evolve_RK4.txt");
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

        V.print("evolve_forward_Euler V= ");
        R.print("evolve_forward_Euler R= ");

        printAndSaveToFile(V, R, "./out/evolve_forward_Euler.txt");
    }

    void printAndSaveToFile(mat V, mat R, string filename)
    {
        int width = 12;
        int prec = 4;

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
