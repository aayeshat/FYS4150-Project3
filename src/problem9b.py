# multidimensional arrays
import numpy as np
import math

# inline plots
import matplotlib.pyplot as plt

# nicer figures
import matplotlib as mpl

t, r1_x, r1_y, r1_z, r2_x, r2_y, r2_z = np.loadtxt(
    "./out/data_R.txt", usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True, skiprows=0
)

plt.plot(r1_x, r1_y)
plt.plot(r2_x, r2_y)
# add axis labels
plt.xlabel("Time [$\mu s$]")
plt.ylabel("Motion in z direction")
# add grid
plt.grid(linestyle="--", linewidth=0.2)
# save plot as pdf
plt.savefig("./out/plot_9b.pdf")



# plt.figure()
# ax = plt.axes(projection="3d")
# plt.tight_layout()
# ax.plot3D(r_x, r_y, r_z, 'blue', label='Trajectory of particle 1 with interactions')
# plt.legend()
# plt.savefig("./out/plot_3d.pdf")
