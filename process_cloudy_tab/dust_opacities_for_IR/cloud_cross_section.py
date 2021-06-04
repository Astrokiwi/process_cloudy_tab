import numpy as np
from scipy import interpolate as interp

Msun = 1.9891e33 # in grams
pc = 3.08567758128e18 # in cm
amu = 1.66054e-24 # in grams
mean_molecular_mass = 1.22


conversion_factor = 1./(2./3./np.pi/mean_molecular_mass/amu*Msun/pc**2)

def full_NH_to_table_density(NH):
    """Convert NH of uniform density cloud to tabulated Msun/pc**2 density.
    NH is hydrogen density through centre of cloud (i.e. through 2*radius)
    Tabulated surface density is mass of cloud/radius**2"""
    return conversion_factor*NH
    

def get_cloud_cross_section_function(file="data/constant_density_table_120319.dat",L=100,dens_bounds=[-9.,8.],temp_bounds=[1.,4.]):
    """ Read density table from file, put it on a regular grid, and return an interpolator
    Interpolater has two inputs: log10 of "surface density", i.e. mass of cloud in solar masses/(radius in pc)**2
    and log10 temperature in K.
    It returns the cross section "factor". Multiply this by the geometric cloud cross-sectional area (i.e. pi * radius**2 )
    to get the true cross section of the cloud.
    """
    indata = np.loadtxt(file,skiprows=1)

    density_measure = np.log10(indata[:,2]/indata[:,1]**2) # log mass/rad**2, mass in Msun, rad in pc
    temps = np.log10(indata[:,0])
    cross_factors = np.log10(indata[:,3]/np.pi/indata[:,1]**2)

    # fix numerical error at large densities - cross_factor can't exceed 1 !!
    cross_factors[cross_factors>0.]=0.

    # log ranges
    dens_range = np.linspace(dens_bounds[0],dens_bounds[1],L)
    temp_range = np.linspace(temp_bounds[0],temp_bounds[1],L)

    p_x,p_y=np.meshgrid(dens_range,temp_range)

    grid_cross = interp.griddata((density_measure,temps),cross_factors,(p_x,p_y))

    return interp.RectBivariateSpline(dens_range,temp_range,grid_cross)
