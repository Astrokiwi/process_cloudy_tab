import numpy as np
import sys

import matplotlib as mpl
mpl.use('Agg')

import pylab as P

insert_dusty = True

# if re-running in same name-space, don't need to reload the data
# this doesn't seem to work now?
if ( not 'd_in' in dir() ):
    print("Loading in table")
    # load in giant file
    #d_in = np.loadtxt("271117.txt",skiprows=1)
    d_in = np.loadtxt("nodust_301117.txt",skiprows=1)
else:
    print("Using table already in memory - hopefully you want to do this!")

d = np.copy(d_in)

# find all values for independent variables (except column density)
denses = np.unique(d[:,1])
intensities = np.unique(d[:,2])
temps = np.unique(d[:,3])

nd = denses.size
nt = temps.size
ni = intensities.size

nc = d.shape[0]//nd//nt//ni
d_offset = nt*ni*nc*np.arange(nd)
i_offset = nt*nc*(np.arange(intensities.size,0,-1)-1) # intensities *decrease* through the array
t_offset = nc*np.arange(temps.size)

if insert_dusty:
    s = d.shape[0]
    d = np.insert(d,5,np.zeros(s),axis=1)
    d = np.insert(d,9,np.zeros(s),axis=1)


#taus = -np.log(d[:,12])
taus = d[:,12]

tau_at_target = np.zeros((ni,nd,nt))

#tau_target = .5

for depth_target in range(16,27):

    for ii in range(ni):
        for id in range(nd):
            for it in range(nt):
                index0 = d_offset[id]+i_offset[ii]+t_offset[it]
                index1 = index0+nc

                tau_at_target[ii,id,it] = taus[np.searchsorted(d[index0:index1,4],10.**depth_target)+index0]

    fig,sps = P.subplots(1,ni,figsize=(2.*ni,4.),dpi=300,sharex=True,sharey=True)

    y,x = np.meshgrid(temps,denses)

    cmin = np.nanmin(tau_at_target)
    cmax = np.nanmax(tau_at_target)

#     cmin = np.percentile(tau_at_target,0)
#     cmax = np.percentile(tau_at_target,90)

    for isp in range(ni):
        sp = sps[isp]
        v=sp.pcolormesh(tau_at_target[isp,:,:].T,cmap='plasma',vmin=cmin,vmax=cmax)
        print(np.min(tau_at_target[isp,:,:]),np.max(tau_at_target[isp,:,:]))
        sp.set_title(r"$I={}$".format(intensities[isp]))
        sp.set_xticks(range(nd))
        sp.set_xlabel(r"$\log n$")
        sp.set_xticklabels(map(str,denses),rotation=90)
        if ( isp==0 ):
            sp.set_ylabel(r"$\log T$")
            sp.set_yticks(range(nt))
            sp.set_yticklabels(map(str,temps))

    P.colorbar(v)
    fig.tight_layout()
    depth_str = "%d"%depth_target
    #P.savefig("../../figures/tau_at_depth"+depth_str+".png")
    P.savefig("../../figures/nodust_tau_at_depth"+depth_str+".png")

    P.close('All')

































