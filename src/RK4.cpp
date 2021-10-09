#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

//Sketch for RK4

using namespace arma;
using namespace std;

void evolve_RK4(double dt){

  int i;
  int t0 = 0;
  int tf; //or put thin is main?
  double dt;
  vec k1, k2, k3, k4, l1, l2, l3, l4;

  vec t = linspace(t0, tf, dt);
  vec a = total_force[i] / particles_.m_; //acc.

  for (int i = 0; i < tf; i++){
        k1 = dt * a[i];
        k2 = dt * a[i]




    }


  }



}
