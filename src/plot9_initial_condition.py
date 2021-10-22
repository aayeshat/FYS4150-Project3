import numpy as np
import math

# inline plots
import matplotlib.pyplot as plt


#  figures
import matplotlib as mpl

mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["legend.fontsize"] = 7.5
plt.rcParams["figure.figsize"] = (7, 5)

interaction = True # True or False

filename=""
outfilename=""

if interaction == True :
  filename= "./out/r_xy_inter_1_2_initial_condition.txt"
  outfilename="./out/plot_9b_inter_initial_condition.pdf"
else:
  filename= "./out/r_xy_nointer_1_2_initial_condition.txt"
  outfilename="./out/plot_9b_nointer_initial_condition.pdf"


t = np.loadtxt(filename, usecols=0, dtype="double")


# with interaction
# reading r values two particles
r_x_1 = np.loadtxt(filename, usecols=1, dtype="double")
r_y_1 = np.loadtxt(filename, usecols=2, dtype="double")
r_z_1 = np.loadtxt(filename, usecols=3, dtype="double")
r_x_2 = np.loadtxt(filename, usecols=4, dtype="double")
r_y_2 = np.loadtxt(filename, usecols=5, dtype="double")
r_z_1 = np.loadtxt(filename, usecols=6, dtype="double")
# Plotting motion of two particles with interactions
plt.scatter(r_x_1, r_y_1, color="blue", label="particle 1")
plt.scatter(r_x_2, r_y_2, color="red", label="particle 2")

# add axis labels and legend
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(linestyle="--", linewidth=0.2)
plt.savefig(outfilename)
