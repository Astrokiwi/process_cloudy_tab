import numpy as np
import sys

import matplotlib as mpl
mpl.use('Agg')

import pylab as P


# if re-running in same name-space, don't need to reload the data
# this doesn't seem to work now?
if ( not 'd_in' in dir() ):
    print("Loading in table")
    # load in giant file
    d_in = np.loadtxt("tables_251017.txt",skiprows=1)
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

nc = d.shape[0]/nd/nt/ni
d_offset = nt*ni*nc*np.arange(nd)
i_offset = nt*nc*(np.arange(intensities.size,0,-1)-1) # intensities *decrease* through the array
t_offset = nc*np.arange(temps.size)

taus = -np.log(d[:,12])

tau_one = np.zeros((ni,nd,nt))

#tau_target = .5

for tau_target in [.5,1.,2.,5.]:

    for ii in range(ni):
        for id in range(nd):
            for it in range(nt):
                index0 = d_offset[id]+i_offset[ii]+t_offset[it]
                index1 = index0+nc

                tau_one[ii,id,it] = d[np.searchsorted(taus[index0:index1],tau_target)+index0,4]
    #             tau_one[ii,id,it] = np.searchsorted(taus[index0:index1],1.)
                #tau_one[ii,id,it] = index0

    tau_one = np.log10(tau_one)

    fig,sps = P.subplots(1,ni,figsize=(2.*ni,4.),dpi=300,sharex=True,sharey=True)

    y,x = np.meshgrid(temps,denses)

    #cmin = np.nanmin(tau_one)
    #cmax = np.nanmax(tau_one)

    #cmin = np.percentile(tau_one,0)
    #cmax = np.percentile(tau_one,90)
    
    cmin = 23.
    cmax = 26.

    for isp in range(ni):
        sp = sps[isp]
        v=sp.pcolormesh(tau_one[isp,:,:].T,cmap='plasma',vmin=cmin,vmax=cmax)
        print(np.min(tau_one[isp,:,:]),np.max(tau_one[isp,:,:]))
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
    tau_str = "%d"%np.floor(tau_target)+"pt"+"%d"%(np.remainder(tau_target,1.)*10.)
    P.savefig("../../figures/tau"+tau_str+".png")

    P.close('All')

































