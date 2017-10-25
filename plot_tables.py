import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')

import pylab as P
#from joblib import Parallel, delayed

# if re-running in same name-space, don't need to reload the data
# this doesn't seem to work now?
if ( not 'd_in' in dir() ):
    print("Loading in table")
    # load in giant file
    d_in = np.loadtxt("tables_161017.txt",skiprows=1)
else:
    print("Using table already in memory - hopefully you want to do this!")


    
print("Working on table")
d = np.copy(d_in) # so we don't repeat the logging, but we can change the logging without having to read in again

# find all values for independent variables (except column density)
denses = np.unique(d[:,1])
intensities = np.unique(d[:,2])
temps = np.unique(d[:,3])

nd = denses.size
nt = temps.size
ni = intensities.size

d[:,4:9] = np.log10(d[:,4:9]) # log most dependent variables, plus column density (other independent variables are already logged) (don't log dust fraction)
#d[:,5:9] = np.log10(d[:,5:9]) # log most dependent variables, BUT NOT column density (other independent variables are already logged) (don't log dust fraction)

d[:,6]-=d[:,1] # convert from erg/s/cm^3 to erg/s/(particle)
d[:,7]-=d[:,1] # convert from erg/s/cm^3 to erg/s/(particle)

# output all different quantities
titles = ["tgrain","heat","cool","prad","dg","kabs","kscat"]
figs = []
sps = []
yranges = [[np.min(d[:,icol+5]),np.max(d[:,icol+5])] for icol in range(len(titles))]
xrange = [np.min(d[:,4]),np.max(d[:,4])]
x_size = xrange[1]-xrange[0]
x_offset = x_size*.1
y_sizes = [yr[1]-yr[0] for yr in yranges]
y_offsets = [ys*.1 for ys in y_sizes]

print(nd,"x",nt,"plots")

print("Calculating order of arrays")
nc = d.shape[0]/nd/nt/ni
d_offset = nt*ni*nc*np.arange(nd)
i_offset = nt*nc*(np.arange(intensities.size,0,-1)-1) # intensities *decrease* through the array
t_offset = nc*np.arange(temps.size)


print("Setting up subplots")
# set up figures
for ifig, title in enumerate(titles):
    print(title)
    fig,sp = P.subplots(1,1,figsize=(48.,12.),dpi=100)
    fig.suptitle(title)
    figs.append(fig)
    sps.append(sp)

def doplot(d,sps,denses,intensities,temps,id,ni,it):
    print(id,it)
    x_left = it*(x_size+x_offset)
    for ii in range(ni):
#         tab_slice = (d[:,1]==denses[id]) & (d[:,2]==intensities[ii]) & (d[:,3]==temps[it])
        index0 = d_offset[id]+i_offset[ii]+t_offset[it]
        index1 = index0+nc
        d_slice = d[index0:index1]
        #d_slice = d_slice[::10,:] # for speed
        
#         for icol,fig in enumerate(figs):
        for icol,sp in enumerate(sps):
            y_bottom = id*(y_sizes[icol]+y_offsets[icol])
            sp.plot(d_slice[:,4]-xrange[0]+x_left,d_slice[:,5+icol]-yranges[icol][0]+y_bottom)

print("Plotting")

#plot all the subplots
for id in range(nd):
    for it in range(nt):
        doplot(d,sps,denses,intensities,temps,id,ni,it)
        for icol,sp in enumerate(sps):
            sp.text(it*(x_size+x_offset),(id+1)*y_sizes[icol]+id*y_offsets[icol],"n{}T{}".format(denses[id],temps[it]))

xticks_left = np.arange(nt)*(x_size+x_offset)
xticks_right = np.arange(nt)*(x_size+x_offset)+x_size

nxticks = 7
nyticks = 5

xticks_each = np.linspace(0,x_size,nxticks)
xticks_all = np.tile(xticks_each,nt)+np.repeat(xticks_left,nxticks)

xtickvalues = d[np.linspace(0,nc-1,nxticks,dtype=int),4]
xticklabels = np.tile(map(lambda x : "%.1f" % x,xtickvalues),nt)

print("Rendering and dumping")
# format & dump everything
for ifig, fig in enumerate(figs):
    sp = sps[ifig]
    sp.set_xlim([0,nt*x_size+nt*x_offset])
    sp.set_ylim([0,nd*y_sizes[ifig]+nd*y_offsets[ifig]])
    sp.set_xticks(xticks_all)
    sp.set_xticklabels(xticklabels,rotation=90)


    yticks_bottom = np.arange(nd)*(y_sizes[ifig]+y_offsets[ifig])
    yticks_each = np.linspace(0,y_sizes[ifig],nyticks)
    yticks_all = np.tile(yticks_each,nd)+np.repeat(yticks_bottom,nyticks)

    ytickvalues = np.linspace(yranges[ifig][0],yranges[ifig][1],nyticks)
    yticklabels = np.tile(map(lambda x : "%.1f"%x,ytickvalues),nd)
    #yticklabels = np.tile(map(lambda x : "{}".format(x),ytickvalues),nd)
    
    sp.set_yticks(yticks_all)
    sp.set_yticklabels(yticklabels)

    fig.tight_layout()
    fig.savefig("../../figures/table_summary_final_"+titles[ifig]+".png")

# to avoid memory getting full of figures
P.close('All')
