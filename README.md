# FYS4150-PROJECT3

Reposistory FYS4130 is created to study the effects for penningtrap for Ca+. Penningtrap is  a device used to store charged particles. This device is known as a Penning trap and utilizes a static configuration of electric and magnetic fields to confine charged particles, so that they can be used for various types of measurements and experiments. In particular, experiments at CERN such as ALPHA, AEgIS and BASE use Penning traps to control the antimatter in their experiments.

In the folder named "Source/Src" there are .hpp and .cpp files for both problem 9 and 10. In folder "Phyton".py files used to make the plots.  The folder "data" have all .txt that are used for generating python plots. In folder 9, the files "Particle.cpp" and "PenningTrap.cpp" which define the main classes. The class "particle" stores the main characteristics of a particle such as its charge, mass, somehow initial position vector and velocity. and the class "penning trap" stores particles and have some functions that evolve the particle through time such as applied magnetic field, potential also ecternal as well as forces from one particle to other. For single particle main9a.cpp which is used to study motion of particle in z direction as time evolves(100microsec). /The files "main_RK4.cpp" and main_RK4.cpp" we evolve only one particle through time using the Foward Euler and the Runge-Kutta method (respectively). With them, we generate the .txt files stored in the data folder files), changing the step-size in the ones used to study the relative error. For main9.cpp having two particles in our trap. Changing different desired parameters we generate different .txt files using runge kutta and forward euler which arestore in the folder data used to  produce the rest of the plots of the folder Plots. 
For problem 10 (folder 10) the aim is to study many particles e.g 100 in a trapp and main10.cpp the code is stimulated for reasonace and varing amplitude for applied potential.
 For building and execute each of the cpp and python is given below
## problem 9

## c++

`g++ src/09/main9a.cpp  src/Particle.cpp src/09/Penningtrap.cpp -larmadillo && ./a.out`

`g++ src/09/main9.cpp src/Particle.cpp src/09/Penningtrap.cpp -larmadillo && ./a.out`

`g++ src/09/main_RK4.cpp src/Particle.cpp src/09/Penningtrap.cpp -larmadillo && ./a.out`

`g++ src/09/main9_FE.cpp src/Particle.cpp src/09/Penningtrap.cpp -larmadillo && ./a.out`

### plot

`python3 src/09/plot9a.py`

`python3 src/09/plot9b.py`

`python3 src/09/3dplot.py`

`python3 src/09/phaseplot.py`

`python3 src/09/relative_error.py`

## problem 10

## c++

`g++ src/10/main10.cpp src/Particle.cpp src/10/Penningtrap.cpp -larmadillo && ./a.out`

### plot

`python3 src/10/plot10.py`
