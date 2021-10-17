# multidimensional arrays
import numpy as np
import math

# inline plots
import matplotlib.pyplot as plt

# nicer figures
import matplotlib as mpl

t, r_x, r_y, r_z = np.loadtxt(
    "./out/data_R.txt", usecols=(0, 1, 2, 3), unpack=True, skiprows=0
)

plt.plot(t, r_z)
# add axis labels
plt.xlabel("Time [$\mu s$]")
plt.ylabel("Motion in z direction")
# add grid
plt.grid(linestyle="--", linewidth=0.2)
# save plot as pdf
plt.savefig("./out/plot_9a.pdf")
