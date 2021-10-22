import numpy as np
import matplotlib.pyplot as plt
from typing import List

#

q = 1   # Charge

B_0 = 96.5 # Magnetic field

m = 40.078 # Mass of Ca+

v_0 = 9.65e8 # applied Potential

d =  10000 # dimensinality 10e4
N = 1000
tf = 100
t = np.linspace(0,tf,N)


# Initial values  for position

x_0 = 1
y_0 = 0
z_0 = 1
# Initial values  for position
v_0 = 0
v_0 =1
v_0 = 0
# Analythical solutions

w_0 = q*(B_0) / m

w_z = 2*q*v_0/ (m * d**2)
w_pos = 0.5 * (w_0 + np.sqrt(w_0**2 - 2 * w_z))
w_neg = = 0.5 * (w_0 - np.sqrt(w_0**2 - 2 * w_z))
A_pos = (v_0 + w_neg * x_0) / (w_neg - w_pos)
A_neg = - (v_0 + w_pos * x_0) / (w_neg - w_pos)

''' x(t) = A_pos*cos(w_pos*t) + A_neg*cos(w_neg*t)
y(t) = A_pos*i*sin(w_pos*t) - A_neg*i*sin(-w_neg*t)
z(t) = z_0*(cos(w_z*t)) + i*sin(w_z*t) '''


def r_exact(t):
    x = A_pos*cos(w_pos*t) + A_neg*cos(w_neg*t)
    y = A_pos*i*sin(w_pos*t) - A_neg*i*sin(-w_neg*t)
    z = z_0*(cos(w_z*t)) + i*sin(w_z*t)
    return r = (x, y, z)

    



#  different stepsizes

textnames =


 
    
    # Analythical solution
        
   
    
  
    # Relative error
    
    
    
    # Graphs
            
    plt.plot()
    
    
# Error convergence rate

ecri = []

for i in  range(1,4):  
    
    arg = np.log10(max_error[i]/max_error[i-1]) / np.log10(h[i]/h[i-1])
    ecri.append( arg  )
    
ecr = sum(ecri) / 4
    
# General characteristics of the graph
    
plt.title("Relative error vs time", fontsize=10)
plt.ylabel(r'$r_{i}$')
plt.xlabel("t (μs)")
plt.legend()
plt.show()
plt.savefig("relative_error.pdf")