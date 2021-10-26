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

interaction = False  # True or False

phase = "z"

filename = ""
filename2 = ""
outfilename = ""

if interaction == True:
    filename = "./out/r_xy_inter_1_2.txt"
    filename2 = "./out/v_xy_inter_1_2.txt"
    outfilename = "./out/phase_" + phase + "_inter.pdf"
else:
    filename = "./out/r_xy_nointer_1_2.txt"
    filename2 = "./out/v_xy_nointer_1_2.txt"
    outfilename = "./out/phase_" + phase + "_nointer.pdf"

# with interaction
# reading r values two particles
r_x_1 = np.loadtxt(filename, usecols=1, dtype="double")
r_y_1 = np.loadtxt(filename, usecols=2, dtype="double")
r_z_1 = np.loadtxt(filename, usecols=3, dtype="double")
r_x_2 = np.loadtxt(filename, usecols=4, dtype="double")
r_y_2 = np.loadtxt(filename, usecols=5, dtype="double")
r_z_2 = np.loadtxt(filename, usecols=6, dtype="double")
## reading v values two particles
v_x_1 = np.loadtxt(filename2, usecols=1, dtype="double")
v_y_1 = np.loadtxt(filename2, usecols=2, dtype="double")
v_z_1 = np.loadtxt(filename2, usecols=3, dtype="double")
v_x_2 = np.loadtxt(filename2, usecols=4, dtype="double")
v_y_2 = np.loadtxt(filename2, usecols=5, dtype="double")
v_z_2 = np.loadtxt(filename2, usecols=6, dtype="double")
# Plotting motion of two particles with interactions

if phase == "x":
    plt.plot(r_x_1, v_x_1, color="blue", label="particle 1")
    plt.plot(r_x_2, v_x_2, color="red", label="particle 2")
elif phase == "y":
    plt.plot(r_y_1, v_y_1, color="blue", label="particle 1")
    plt.plot(r_y_2, v_y_2, color="red", label="particle 2")
elif phase == "z":
    plt.plot(r_z_1, v_z_1, color="blue", label="particle 1")
    plt.plot(r_z_2, v_z_2, color="red", label="particle 2")

plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
    # add axis labels and legend
plt.xlabel("" + phase + " (μm) ")
plt.ylabel("v_" + phase + " (μm/s) ")

plt.legend()
plt.grid(linestyle="--", linewidth=0.2)
plt.savefig(outfilename)
