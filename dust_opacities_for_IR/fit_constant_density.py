import numpy as np
from scipy import interpolate as interp
import matplotlib.pyplot as plt

indata = np.loadtxt("data/constant_density_table_120319.dat",skiprows=1)

density_measure = np.log10(indata[:,2]/indata[:,1]**2) # log mass/rad**2, mass in Msun, rad in pc
temps = np.log10(indata[:,0])
cross_factors = np.log10(indata[:,3]/np.pi/indata[:,1]**2)

# fix numerical error at large densities - cross_factor can't exceed 1 !!
cross_factors[cross_factors>0.]=0.

# log ranges
Lgrid=100
dens_range = np.linspace(-9.,8.,Lgrid)
temp_range = np.linspace(1.,4.,Lgrid)

p_x,p_y=np.meshgrid(dens_range,temp_range)

grid_cross = interp.griddata((density_measure,temps),cross_factors,(p_x,p_y))

cross_spline = interp.RectBivariateSpline(dens_range,temp_range,grid_cross)

z_range = [np.min(grid_cross),np.max(grid_cross)]

plt.figure()
plt.pcolormesh(p_x,p_y,grid_cross)
plt.ylabel(r'$\log T$ (K)')
plt.xlabel(r'$\log\Sigma$ (M$_\odot$pc$^{-2}$)')
plt.colorbar(label="$\log \sigma/A$")
plt.contour(p_x,p_y,grid_cross,levels=np.linspace(z_range[0],z_range[-1],20),colors='black',alpha=.5)
plt.savefig("figs/fit_test.png")

