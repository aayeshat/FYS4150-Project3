# multidimensional arrays
import numpy as np
import math

# inline plots
import matplotlib.pyplot as plt


# nicer figures
import matplotlib as mpl
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 7.5
plt.rcParams["figure.figsize"] = (7,5)
#reading time t
t = np.loadtxt("data.txt", usecols=0, dtype='double')

#reading r values
r_z = np.loadtxt("data.txt", usecols=1, dtype='double')

#reading v values
#v_x = np.loadtxt("data.txt", usecols=2, dtype='double')
# Plotting motion in the z direction as a function of time
plt.plot(t, r_z)


# add axis labels
plt.xlabel("Time [$\mu s$]")
plt.ylabel("Motion in z direction")

# add grid
plt.grid(linestyle = '--', linewidth = 0.2)

# save plot as pdf
plt.savefig("plot_9a.pdf")