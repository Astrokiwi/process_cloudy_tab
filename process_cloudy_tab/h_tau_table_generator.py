#!/usr/bin python3


import numpy as np
import scipy.integrate as integrate

def veckernel(x_in):
    x = np.abs(x_in)
    kern = np.zeros(x.shape)
    
    centre_kernel = (x<.5)
    kern[centre_kernel] = 1.-6.*x[centre_kernel]**2+6.*x[centre_kernel]**3
    
    outer_kernel = (x>=.5) & (x<=1.)
    kern[outer_kernel] = 2.*(1.-x[outer_kernel])**3
    
    kern*=(8./np.pi)
    return kern


def kernel(x_in):
    x = np.abs(x_in)
    if x<0.5:
        kern = 1.-6.*x**2+6.*x**3
    elif x<=1.:
        kern = 2.*(1.-x)**3
    else:
        kern = 0.
    kern*=(8./np.pi)
    return kern

def rkernel(r,z):
    x = np.sqrt(r**2+z**2)
    return kernel(x)

nr = 1001
nz = 2001

rcoords = np.linspace(0.,1.,nr)
zcoords = np.linspace(-1.,1.,nz)

rz_mesh = np.meshgrid(rcoords,zcoords)

rad2d_p = rz_mesh[0]
z_p = rz_mesh[1]

rad3d_p = np.sqrt(rad2d_p**2+z_p**2)

kern_p = veckernel(rad3d_p)

centre_max_surface_density = integrate.quad(lambda z: kernel(z),-1.,1.)[0]

cum_surf = integrate.cumtrapz(kern_p[:,:],zcoords,axis=0)
cum_surf/=centre_max_surface_density

nsamples = 101
surf_samples = np.linspace(0.,1.,nsamples)

surf_locs = np.apply_along_axis(lambda x: np.searchsorted(x,surf_samples),axis=0,arr=cum_surf)
zcent_at_sample = zcoords[surf_locs[:,0]]
zcent_at_sample = (zcent_at_sample[1:] + zcent_at_sample[:-1])/2.

surf_locs[surf_locs==cum_surf.shape[0]]=cum_surf.shape[0]-1

cum_surf_at_samples = cum_surf[surf_locs.flatten(),np.tile(np.arange(surf_locs.shape[1]),surf_locs.shape[0])].reshape(surf_locs.shape[0],surf_locs.shape[1])

cum_mass_at_samples = np.sum(cum_surf_at_samples,axis=1)

cum_mass_at_samples/=cum_mass_at_samples[-1]

mass_per_slice = np.diff(cum_mass_at_samples)
sample_centres = (surf_samples[:-1]+surf_samples[1:])/2.

#np.savetxt("data/h_depth_table.dat",np.array([surf_samples[:-1],mass_per_slice,zcent_at_sample[:-1]]).T)
np.savetxt("../data/h_depth_table.dat",np.array([sample_centres,mass_per_slice,zcent_at_sample]).T)


