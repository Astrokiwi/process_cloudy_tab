import numpy as np
import time
import sys

#taumode=True means use tau as the independent variable
#taumode=False means use NH column density instead
taumode = True

#nodustmode = this is dust-free gas for the secondary table
nodustmode = False

#highdensemode = this is the third table, high density but only cold gas
highdensemode = True

# if re-running in same name-space, don't need to reload the data
if ( not 'd_in' in dir() ):
    print("Loading in table")
    # load in giant file
    #d_in = np.loadtxt("tables_271117.txt",skiprows=1)
#     d_in = np.loadtxt("nodust_301117.txt",skiprows=1)
    d_in = np.loadtxt("highden_260118.txt",skiprows=1)
else:
    print("Using table already in memory - hopefully you want to do this!")

d = np.copy(d_in)

if nodustmode:
    # add in blank data
    s = d.shape[0]
    d = np.insert(d,5,np.zeros(s),axis=1)
    d = np.insert(d,9,np.zeros(s),axis=1)

#convert from exp(-tau) to tau
#d[:,12] = -np.log(d[:,12]) # should already be done now # no, it's not logged, but that's fine.

denses = np.unique(d[:,1])
intensities = np.unique(d[:,2])
temps = np.unique(d[:,3])

if taumode:
    column_in = 10.**np.linspace(-2.,.7,50)
else:
    column_in = 10.**np.linspace(18.477,25.,50)

#taus = np.linspace(0.,5.,50)

nd = denses.size
nt = temps.size
ni = intensities.size

tau_suffix = "tau" if taumode else "coldens"
dust_suffix = "nodust" if nodustmode else ""
if nodustmode and highdensemode:
    raise Exception("Can't have both nodustmode and highdensemode!")
dense_suffix = "dense" if highdensemode else ""
 
suffixes = tau_suffix+dust_suffix+dense_suffix
f = open("shrunk_table_labels_"+time.strftime("%d%m%y")+suffixes+".dat",'w')

for ar in [denses,temps,intensities,column_in]:
    #outd = np.hstack([outd,ar.size,ar])
    f.write(str(ar.size)+"\n")
    for v in ar:
        f.write(str(v)+"\n")
    

f.close()

if taumode:
    outp_cols_interp = [6,7,5,8,9,10,11,4]
    icol_in = 12
    outp_dolog = [True,True,False,True,False,False,False,False]
else:
    outp_cols_interp = [6,7,5,8,9,10,11,12]
    icol_in = 4
    outp_dolog = [True,True,False,True,False,False,False,False]
noutp = len(outp_cols_interp)

alloutp = [np.empty((nd,nt,ni),dtype=object) for x in range(noutp)]

nc = d.shape[0]//nd//nt//ni
d_offset = nt*ni*nc*np.arange(nd)
i_offset = nt*nc*(np.arange(intensities.size,0,-1)-1) # intensities *decrease* through the array
t_offset = nc*np.arange(temps.size)


for id in range(nd):
    for it in range(nt):
        print(id,it)
        for ii in range(ni):
            index0 = d_offset[id]+i_offset[ii]+t_offset[it]
            index1 = index0+nc
            d_slice = d[index0:index1]

            for i,icol in enumerate(outp_cols_interp):
                alloutp[i][id,it,ii] = np.interp(column_in,d[index0:index1,icol_in],d[index0:index1,icol])


alloutdata = np.empty((noutp,0))

for id in range(nd):
    for it in range(nt):
        print(id,it)
        for ii in range(ni):
            outdata = []
            for iout,outp in enumerate(alloutp):
                if ( outp_dolog[iout] ):
                    outdata.append(np.log10(outp[id,it,ii]))
                else:
                    outdata.append(outp[id,it,ii])
            outdata = np.array(outdata)
            alloutdata=np.hstack([alloutdata,outdata])

alloutdata = np.array(alloutdata)
np.savetxt("shrunk_table_"+time.strftime("%d%m%y")+suffixes+".dat",alloutdata.T)
