import numpy as np
import math

# inline plots
import matplotlib.pyplot as plt


# nicer figures
import matplotlib as mpl

mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["legend.fontsize"] = 7.5
plt.rcParams["figure.figsize"] = (7, 5)

interaction = True # True or False

filename=""

filename2=""
outfilename=""

if interaction == True :
  filename= "./out/r_xy_inter_1_2.txt" 
  filename2 = "./out/v_xy_inter_1_2.txt"
  outfilename="./out/phase_z_inter.pdf"
else:
  filename= "./out/v_xy_nointer_1_2.txt" 
  filename2  = "./out/v_xy_nointer_1_2.txt"
  outfilename = "./out/phase_z_nointer.pdf"


t = np.loadtxt(filename, usecols=0,  dtype="double")
t2 = np.loadtxt(filename2, usecols=0,  dtype="double")

# with interaction
# reading r values two particles
r_x_1 = np.loadtxt(filename, usecols=1, dtype="double")
r_y_1 = np.loadtxt(filename, usecols=2, dtype="double")
r_z_1 = np.loadtxt(filename, usecols=3, dtype="double")
r_x_2 = np.loadtxt(filename, usecols=4, dtype="double")
r_y_2 = np.loadtxt(filename, usecols=5, dtype="double")
r_z_2 = np.loadtxt(filename, usecols=6, dtype="double")
v_x_1 = np.loadtxt(filename2, usecols=1, dtype="double")
v_y_1 = np.loadtxt(filename2, usecols=2, dtype="double")
v_z_1 = np.loadtxt(filename2, usecols=3, dtype="double")
v_x_2 = np.loadtxt(filename2, usecols=4, dtype="double")
v_y_2 = np.loadtxt(filename2, usecols=5, dtype="double")
v_z_2 = np.loadtxt(filename2, usecols=6, dtype="double")
# Plotting motion of two particles with interactions
plt.plot(r_z_1,v_z_1 , color="blue", label="particle 1")
plt.plot(r_z_2, v_z_2, color="red", label="particle 2")

# add axis labels and legend
plt.xlabel("z")
plt.ylabel("v_z")
plt.legend()
plt.grid(linestyle="--", linewidth=0.2)
plt.savefig(outfilename)
