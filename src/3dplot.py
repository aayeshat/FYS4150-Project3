from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection='3d')
interaction = True # True or False

filename=""
outfilename=""


if interaction == True :
  filename= "./out/r_xy_inter_1_2.txt" 
  outfilename="./out/3d_inter.pdf"
else:
  filename= "./out/r_xy_nointer_1_2.txt" 
  outfilename="./out/3d_nointer.pdf"


t = np.loadtxt(filename, usecols=0, dtype="double")


# with interaction
# reading r values two particles
r_x_1 = np.loadtxt(filename, usecols=1, dtype="double")
r_y_1 = np.loadtxt(filename, usecols=2, dtype="double")
r_z_1 = np.loadtxt(filename, usecols=3, dtype="double")
r_x_2 = np.loadtxt(filename, usecols=4, dtype="double")
r_y_2 = np.loadtxt(filename, usecols=5, dtype="double")
r_z_2 = np.loadtxt(filename, usecols=6, dtype="double")
# Plotting motion of two particles with interactions
plt.plot(r_x_1, r_y_1,r_z_1,  color="blue", label="particle 1")
plt.plot(r_x_2, r_y_2,r_z_2, color="red", label="particle 2")

# add axis labels and legend
plt.xlabel("x")
plt.ylabel("y")
##plt.zlabel("z")
plt.legend()
plt.grid(linestyle="--", linewidth=0.2)
plt.savefig(outfilename)

