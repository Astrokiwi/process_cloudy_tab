#!/usr/bin/python3

import numpy as np
import scipy.integrate as integrate
import time

from sys import path
path.append("../src/")
import tab_interp

print("Loading table (coldens)")
chTab = tab_interp.CoolHeatTab(("shrunk_table_labels_171117coldens.dat"),("shrunk_table_171117coldens.dat"))
interpTabVec = np.vectorize(chTab.interpTab)

print("Loading table (tau)")
tauChTab = tab_interp.CoolHeatTab(("shrunk_table_labels_171117tau.dat"),("shrunk_table_171117tau.dat"))
interpTauTabVec = np.vectorize(tauChTab.interpTab)

def get_table_labels(fname):
    table_labels = np.loadtxt(fname)

    agn_ndense = int(table_labels[0])
    agn_dense_vals = table_labels[1:1+agn_ndense]

    agn_ntemp = int(table_labels[1+agn_ndense])
    agn_temp_vals = table_labels[2+agn_ndense:2+agn_ndense+agn_ntemp]

    agn_nintensity = int(table_labels[2+agn_ndense+agn_ntemp])
    agn_intensity_vals = table_labels[3+agn_ndense+agn_ntemp:3+agn_ndense+agn_ntemp+agn_nintensity]

    agn_ncolumn_in = int(table_labels[3+agn_ndense+agn_ntemp+agn_nintensity])
    agn_column_in_vals = table_labels[4+agn_ndense+agn_ntemp+agn_nintensity:]
    
    return agn_ndense,agn_dense_vals,agn_ntemp,agn_temp_vals,agn_nintensity,agn_intensity_vals,agn_ncolumn_in,agn_column_in_vals

agn_ndense,agn_dense_vals,agn_ntemp,agn_temp_vals,agn_nintensity,agn_intensity_vals,agn_ncolumn_in,agn_column_in_vals = get_table_labels("shrunk_table_labels_171117coldens.dat")

*x,agn_ntau_in,agn_tau_vals = get_table_labels("shrunk_table_labels_171117tau.dat")

print("Tables loaded")

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

get_struct_attribute = np.vectorize(lambda x,y:x.__getattribute__(y))

h_tau_table = np.loadtxt("../data/h_depth_table.dat")
tab_surfs = h_tau_table[:,0]
tab_dm = h_tau_table[:,1]
tab_dz = h_tau_table[:,2]

# for test
#tab_surfs = np.array([0.,0.])
#tab_dm = np.array([.5,.5])

#attributes_to_mean = ['dHeat','dCool','dg','dustT','arad','opac_abs','opac_scat','column_out'] # for production - i.e. output by optical depth, for raytracing better
attributes_to_mean = ['dHeat','dCool','dustT','arad','dg','opac_abs','opac_scat','column_in'] # for test - i.e. same as input, output by column density

centre_max_surface_density = integrate.quad(lambda z: kernel(z),-1.,1.)[0]

# h=nH**(-1./3.)/h_dense_factor
# h comes out in pc
h_dense_factor = 0.3404

# conversion factor from solar masses/pc**2 into N_H in cm**-2
molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
solar_mass_cgs = 1.989e33
pc_cgs = 3.086e18

mass_density_to_number_density = solar_mass_cgs/pc_cgs**2/(molecular_mass*proton_mass_cgs)

mass_p = 0.1 # solar masses

#nH_p = 1.e4
#nH_p = 10.**np.linspace(1.,7.,7)
nH_p = 10.**np.repeat(agn_dense_vals,agn_ntemp*agn_nintensity*agn_ncolumn_in)
TK_p = 10.**np.tile(np.repeat(agn_temp_vals,agn_nintensity*agn_ncolumn_in),agn_ndense)
flux_p = 10.**np.tile(np.repeat(agn_intensity_vals,agn_ncolumn_in),agn_ndense*agn_ntemp)
surf0_p = np.tile(agn_column_in_vals,agn_nintensity*agn_ndense*agn_ntemp)

surf0_p = np.broadcast_to(surf0_p,(tab_surfs.size,surf0_p.size)).T

# TK_p = np.full_like(nH_p,500.)
# flux_p = np.full_like(nH_p,10.**5)

h_p = nH_p**(-1./3.)*h_dense_factor
surface_density = centre_max_surface_density*mass_p/h_p**2
surface_number_density = surface_density * mass_density_to_number_density

surfs_interior = np.outer(surface_number_density,tab_surfs)+surf0_p
#TK_p = np.full_like(surfs_interior,500.)
#flux_p = np.full_like(surfs_interior,10.**5)
#nH_p = np.repeat(nH_p.reshape(7,1),100,1)

nH_p = np.broadcast_to(nH_p,(tab_surfs.size,nH_p.size)).T
TK_p = np.broadcast_to(TK_p,(tab_surfs.size,TK_p.size)).T
flux_p = np.broadcast_to(flux_p,(tab_surfs.size,flux_p.size)).T


tabStructs = interpTabVec(nH_p.astype(np.float64),TK_p.astype(np.float64),flux_p.astype(np.float64),surfs_interior.astype(np.float64))


attrib_out = []
for attribute in attributes_to_mean:
#for attribute in ['dustT']:
    if attribute=='column_in':
        tau_indices = get_struct_attribute(tabStructs,'ic')
        vals = agn_column_in_vals[tau_indices]
        print(vals.shape)
    else:
        vals = get_struct_attribute(tabStructs,attribute)
        print(vals.shape)
    vals = vals * tab_dm
    vals = np.sum(vals,1)
    attrib_out.append(vals)
    
attrib_out = np.array(attrib_out)

np.savetxt("shrunk_table_"+time.strftime("%d%m%y")+"_hsmooth_tau.dat",attrib_out.T)
