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

interaction = False  # True or False

filename1 = ""
filename2 = ""
filename3 = ""
outfilename = ""


if interaction == True:
    filename1 = "./out/data/10_inter__0.100000.txt"
    filename2 = "./out/data/10_inter__0.400000.txt"
    filename3 = "./out/data/10_inter__0.700000.txt"
    outfilename = "./out/plot/plot_10_inter.pdf"
else:
    filename1 = "./out/data/10_nointer__0.100000.txt"
    filename2 = "./out/data/10_nointer__0.400000.txt"
    filename3 = "./out/data/10_nointer__0.700000.txt"
    outfilename = "./out/plot/plot_10_nointer.pdf"


omega = np.loadtxt(filename1, usecols=(0), delimiter="  ", dtype="double")

n_f1 = np.loadtxt(filename1, usecols=(1), delimiter="  ", dtype="double")
n_f4 = np.loadtxt(filename2, usecols=(1), delimiter="  ", dtype="double")
n_f7 = np.loadtxt(filename3, usecols=(1), delimiter="  ", dtype="double")

plt.plot(omega, n_f1, label="f=0.1", color="c")
plt.plot(omega, n_f4, label="f=0.4", color="y")
plt.plot(omega, n_f7, label="f=0.7", color="b")
plt.legend()

plt.xlabel("$\omega_v [MHz]$")
plt.ylabel("Fraction of particle ")

plt.grid(linestyle="--", linewidth=0.2)

plt.savefig(outfilename)
