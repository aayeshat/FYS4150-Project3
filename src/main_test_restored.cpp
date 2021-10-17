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
    double t0 = 0;

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
  vec coulomb_interaction(int i, int j)
    {
        if (i != j)
        {
            return force_particle(i, j);
        }
        else
        {
            return vec(3).fill(0);
        }
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

       for (int j = 0; j < n + 1; j++)
      {
           if (i != j)
            {
                F += coulomb_interaction(i, j);
            }
           
        }
        return F;
    }

    vec total_force(int i)
    {
        //vec F = total_force_particles(i) + total_force_external(i);
        vec F = total_force_external(i);

        return F;
    }

    void evolve_RK4(double dt)
    {
        mat R = mat(3, n).fill(0);
        mat V = mat(3, n).fill(0);

        for (int i = 0; i < n + 1; i++)
        {
            Particle p = particles[i];

            vec F = total_force(i);
            vec acceleration = F / p.m;

            vec v_0 = p.v; // particle initial velocity v0
            vec r_0 = p.r; // particle initial position

            // K1 velocity and position
            vec k1_v = acceleration * dt;

            vec k1_r = v_0 * dt;

            // K2
            vec k2_v = dt * (t0 + dt / 2, v_0 + k1_v / 2);

            vec k2_r = dt * (t0 + dt / 2, r_0 + k1_r / 2);

            double t_1 = t0;
            vec v_1 = v_0;
            vec r_1 = r_0;

            // K3
            vec k3_v = dt * (t_1 + dt / 2, v_1 + k2_v / 2);

            vec k3_r = dt * (t_1 + dt / 2, r_1 + k1_r / 2);

            double t_2 = t_1;
            vec v_2 = v_1;
            vec r_2 = r_1;

            // K4
            vec k4_v = dt * (t_2 + dt / 2, v_2 + k3_v);

            vec k4_r = dt * (t_2 + dt / 2, r_2 + k3_r);

            V.col(i) = v_0 + (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6;
            R.col(i) = r_0 + (k1_r + 2 * k2_r + 2 * k3_r + k4_r) / 6;
        }

        V.print("evolve_RK4 V= ");
        R.print("evolve_RK4 R= ");

        printAndSaveToFile(V, R, "./out/evolve_RK4.txt");
    }

    void evolve_forward_Euler(double dt)
    {
        mat R = mat(3, n).fill(0);
        mat V = mat(3, n).fill(0);

        for (int i = 0; i < n ; i++)
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

int main()
{

    double B0 = 9.65e1;
    double V0 = 9.65e8;
    double d = 10e4;
    int n = 2; //number of particles

    int t = 100;
    int N = 1000;
    
    double dt = t * (1. / N);
    //vec t = linspace(t0, tf, dt);

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

   /*  cout << endl
         << endl
         << "evolve_forward_Euler "
         << endl
         << endl;
    pt.evolve_forward_Euler(.1); */

    return 0;
}
