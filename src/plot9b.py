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

t = np.loadtxt("./out/xy_inter_1_2.txt", usecols=0, dtype="double")


# with interaction
# reading r values two particles
r_x_1 = np.loadtxt("./out/xy_inter_1_2.txt", usecols=1, dtype="double")
r_y_1 = np.loadtxt("./out/xy_inter_1_2.txt", usecols=2, dtype="double")

r_x_2 = np.loadtxt("./out/xy_inter_1_2.txt", usecols=4, dtype="double")
r_y_2 = np.loadtxt("./out/xy_inter_1_2.txt", usecols=5, dtype="double")

# Plotting motion of two particles with interactions
plt.plot(r_x_1, r_y_1, color="blue", label="particle 1")
plt.plot(r_x_2, r_y_2, color="red", label="particle 2")

# add axis labels and legend
plt.xlabel("x")
plt.ylabel("y")
plt.legend()

# add grid
plt.grid(linestyle="--", linewidth=0.2)

# save plot as pdf
plt.savefig("./out/plot_9b.pdf")
