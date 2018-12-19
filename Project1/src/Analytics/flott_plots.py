import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from data_analytics import read_parameters
sns.set_style("darkgrid")

"""
This is for plotting the onebody density. This is not stream lined, so
one has to make a density and density_non file (hardcoding the name into the
cpp file simulation) if one wants to compare them.

This will also calculate the mean r as seen in the rapport.

The reason N is multiplied is that we have normalized rho(r) to 1
in out calculation in simulation.cpp
"""

parameters = read_parameters("../output/")
N = parameters["N"]

r = np.fromfile("../output/r_positions.bin",sep=" ")
rho = np.fromfile("../output/density.bin",sep=" ")
#rho_non = np.fromfile("../output/density_non.bin",sep=" ")
volume_factor = np.fromfile("../output/volume.bin",sep=" ")


rs = np.linspace(r[0],r[1],int(r[2]),endpoint=True)
step = r[3]

mean = np.sum(rho[:len(rs)]*volume_factor[:len(rs)]*rs)
print("Mean r: ",mean)

plt.plot(rs,N*rho[:len(rs)],label="Interacting")
#plt.plot(rs,N*rho_non[:len(rs)],label="Noninteracting")
plt.title(r"Onebody density for N=%g"%N,fontsize=15)
plt.xlabel("$r$",fontsize=20)
plt.ylabel(r"$\rho(r)$",fontsize=20)
plt.xlim(0,4)
plt.legend()
plt.show()
