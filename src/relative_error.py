import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

q = 1   # Charge
B0 = 96.5 # Magnetic field
m = 40.078 # Mass of Ca+
V0_d2 = 9.65 # Applied potential over d squared
N = 100000 # Steps
tf = 100 #Final time
d = 1e4
t = np.linspace(0, tf, N)

# Initial values
x_0 = 10
y_0 = 0
z_0 = 10
v_0 = 10 #velocity in y direction

# Variables
w_0 = q*(B0) / m
w_z = np.sqrt(2 * q * V0_d2 / m)  #omega z squared
w_pos = 0.5 * (w_0 + np.sqrt(w_0**2 - 2 * w_z)) # omega+
w_neg = 0.5 * (w_0 - np.sqrt(w_0**2 - 2 * w_z)) # omega-
A_pos = (v_0 + w_neg * x_0) / (w_neg - w_pos) # A+
A_neg = - (v_0 + w_pos * x_0) / (w_neg - w_pos) # A-

#Analytical solutions


filename_rk = ['../out/r_1(0.1).txt', '../out/r_1(0.05).txt', '../out/r_1(0.01).txt', '../out/r_1(0.005).txt', '../out/r_1(0.001).txt']
h_name = ['0.1', '0.05', '0.01', '0.005','0.001']
filename_fe = ['../out/r_1_euler(0.1).txt','../out/r_1_euler(0.05).txt','../out/r_1_euler(0.01).txt','../out/r_1_euler(0.005).txt','../out/r_1_euler(0.001).txt']

#Loop for plotting relative error for 5 different stepsizes h

max_error_rk = []
max_error_fe = []

for i in range(5):
    x_rk, y_rk, z_rk = np.loadtxt(filename_rk[i], usecols = (1, 2, 3), unpack = True)
    x_fe, y_fe, z_fe = np.loadtxt(filename_fe[i], usecols = (1, 2, 3), unpack = True)

    t = np.linspace(0,tf, len(x_rk))

    x_exc = A_pos * np.cos(w_pos * t) + A_neg * np.cos(w_neg * t)
    y_exc = - A_pos * np.sin(w_pos * t) - A_neg * np.sin(w_neg * t)
    z_exc = z_0 * (np.cos(w_z * t)) + np.sin(w_z * t)

    r_exc = np.sqrt(x_exc**2 + y_exc**2 + z_exc**2)
    r_rk = np.sqrt(x_rk**2 + y_rk**2 + z_rk**2)
    r_fe = np.sqrt(x_fe**2 + y_fe**2 + z_fe**2)

    rel_err_rk = np.absolute( r_exc - r_rk ) / (r_exc)
    rel_err_fe = np.absolute( r_exc - r_fe ) / (r_exc)

    max_rk = max(rel_err_rk)
    max_fe = max(rel_err_fe)

    max_error_rk.append(max_rk)
    max_error_fe.append(max_fe)


    # Plot RK solutions for different stepsizes

    # plt.plot(t, rel_err_rk, label = "h = " + h_name[i])
    plt.ylabel('Relative error $\epsilon_i$')
    plt.xlabel('Time / $\mu$s')
    # plt.savefig('../out/relerr_rk.pdf')

    # Plot FE solutions for different stepsizes

    plt.plot(t, rel_err_fe, label = "h = " + h_name[i])
    plt.legend()
    plt.savefig('../out/relerr_fe.pdf')


#Plot single numerical solution

# x_rk, y_rk, z_rk = np.loadtxt('../out/r_1_euler(0.001).txt', usecols = (1, 2, 3), unpack = True)
# plt.plot(t,x_rk, label = 'x')
# plt.plot(t,y_rk, label = 'y')
# plt.plot(t,z_rk, label = 'z')
# plt.ylabel('x, y, z / $\mu$m')
# plt.xlabel('time / $\mu$s')
# plt.legend()
# plt.savefig('../out/numerical_solutions.pdf')
# plt.show()

h = [0.1, 0.05, 0.01, 0.005, 0.001]

r_err_rk = []
r_err_fe = []


for k in range(1,5):
    h_den = (np.log10(h[k] / h[k-1]))

    rk = 0.25 * np.log10( max_error_rk[k] / max_error_rk[k-1] ) / h_den
    r_err_rk.append(rk)

    fe = 0.25 * np.log10( max_error_fe[k] / max_error_fe[k-1] ) / h_den
    r_err_fe.append(fe)

# Plot estimation of error convergence

# plt.plot(r_err_rk, marker = 'o', label = 'RK')
# plt.plot(r_err_fe, marker = 'o', label = 'FE')
# plt.xlabel('k - 2')
# plt.ylabel('$r_{err}$')
# plt.xticks([0,1,2,3,])
# plt.legend()
#plt.savefig('../out/error_conv.pdf')
# plt.show()
