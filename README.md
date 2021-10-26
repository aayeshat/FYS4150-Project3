# FYS4150-PROJECT3

Reposistory FYS4130 is created to study the effects of ODE for penningtrap. Penningtrap is  a device used to store charged particles. This device is known as a Penning trap and utilizes a static configuration of electric and magnetic fields to confine charged particles, so that they can be used for various types of measurements and experiments. In particular, experiments at CERN such as ALPHA, AEgIS and BASE use Penning traps to control the antimatter in their experiments.

In the folder named "Source/Src" there are .hpp and .cpp files. In folder "Phyton".py files used to make the plots.  The folder "data" have all .txt that are used for generating python plots.

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
