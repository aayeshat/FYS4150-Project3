import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["legend.fontsize"] = 7.5
plt.rcParams["figure.figsize"] = (7, 5)

t = np.loadtxt("./out/r_1.txt", usecols=(0), delimiter='  ', dtype="double")
r_z = np.loadtxt("./out/r_1.txt", usecols=(3), delimiter='  ', dtype="double")

plt.plot(t, r_z)

plt.xlabel("Time [$\mu s$]")
plt.ylabel("Motion of single particle in z direction")

plt.grid(linestyle = '--', linewidth = 0.2)

plt.savefig("../out/plot_9a.pdf")
