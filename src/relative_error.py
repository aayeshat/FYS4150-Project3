import numpy as np
import matplotlib.pyplot as plt
from glob import glob

q = 1   # Charge
B0 = 96.5 # Magnetic field
m = 40.078 # Mass of Ca+
V0_d2 = 9.65 # Applied potential over d squared
N = 100000 # Steps
tf = 100 #Final time
t = np.linspace(0,tf,N) #time

# Initial values
x_0 = 1
y_0 = 0
z_0 = 1
v_0 = 1 #velocity in y direction

# Variables
w_0 = q*(B0) / m
w_z = 2 * q * V0_d2 / m  #omega z squared
w_pos = 0.5 * (w_0 + np.sqrt(w_0**2 - 2 * w_z)) # omega+
w_neg = 0.5 * (w_0 - np.sqrt(w_0**2 - 2 * w_z)) # omega-
A_pos = (v_0 + w_neg * x_0) / (w_neg - w_pos) # A+
A_neg = - (v_0 + w_pos * x_0) / (w_neg - w_pos) # A-

#Analytical solutions
x_exc = A_pos * np.cos(w_pos * t) + A_neg * np.cos(w_neg * t)
y_exc = A_pos * np.sin(w_pos * t) - A_neg * np.sin(w_neg * t)
z_exc = z_0 * (np.cos(w_z * t)) + np.sin(w_z * t)

# h_names = [ "0.01" , "0.05" , "0.1" , "0.5" , "1"]
#Filenames for numerical output from RK
# filenames = ['', '', '', '', '']
# filename = ''

#Loop for plotting relative error for 5 different stepsizes h
#for i in range(4):

    #Load x y z columns from output files

#Load x y z output from RK
fname_rk = '../out/r_xy_nointer_1_2.txt'
x, y, z = np.loadtxt(fname_rk, usecols = (1,2,3), unpack = True)

r_exc = np.sqrt(x_exc**2 + y_exc**2 + z_exc**2)*100 #multiply with 100 for numbers on same order
r_rk = np.sqrt(x**2 + y**2 + z**2)

#Relative error RK
rel_err = ( r_exc - r_rk ) / r_exc

#Load x y z output from FE
fname_euler = '../out/r_xy_nointer_1_2_euler.txt'
x_fe, y_fe, z_fe = np.loadtxt(fname_euler, usecols = (1,2,3), unpack = True)

#Relative error FE
r_fe = np.sqrt(x_fe**2 + y_fe**2 + z_fe**2)
rel_err_fe = ( r_exc - r_fe ) / r_exc

# Plotting relative error
plt.plot(t, rel_err, label = 'RK')
plt.plot(t, rel_err_fe, label = 'FE')
plt.legend()
plt.ylabel('Relative error')
plt.xlabel('t')
plt.savefig('../out/rel_err.pdf')
plt.show()


#Plotting analytical solutions
plt.plot(t, x_exc, label = 'x')
plt.plot(t, y_exc, label = 'y')
plt.plot(t, z_exc, label = 'z')
plt.ylabel('x,y,z')
plt.xlabel('t')
plt.legend()
plt.savefig('../out/analytical_solutions.pdf')


#Plotting numerical solutions
plt.plot(t, x, label = 'x')
plt.plot(t, y, label = 'y')
plt.plot(t, z, label = 'z')
plt.ylabel('x,y,z')
plt.xlabel('t')
plt.legend()
plt.savefig('../out/numerical_solutions.pdf')
