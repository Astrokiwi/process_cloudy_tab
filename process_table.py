import numpy as np
import time
import sys

#d = np.loadtxt("cooling_total.txt",skiprows=1)
#d = np.loadtxt("dust2.txt",skiprows=1)
#d = np.loadtxt("a_rad.txt",skiprows=1)
#d = np.loadtxt("dg.txt",skiprows=1)
#d = np.loadtxt("opac.txt",skiprows=1)
#d = np.loadtxt("cooling_final.txt",skiprows=1)
d = np.loadtxt("tables_021117.txt",skiprows=1)

# intensities = np.linspace(3.748188,7.748188,num=9)
# # temps = np.linspace(1.,8.,29)
# temps = np.array(
#        (1.  ,  1.25,  1.5 ,  1.75,  2.  ,  2.25,  2.5 ,  2.75,  3.  ,
#         3.25,  3.5 ,  3.75,  4.  ,  4.25,  4.5 ,  4.75,  5.  ,  5.25,
#         5.5 ,  5.75,  6.  ,  7.  ,  8.  )
#         )
# denses = np.linspace(2.,9.,8)
# 
# surfs = 10**np.linspace(16.,25.,50)
# 
# nd = denses.size
# nt = temps.size
# ni = intensities.size
# nn = surfs.size

# print(np.unique(d[:,1]),np.unique(d[:,2]),np.unique(d[:,3]))
# print(denses,intensities,temps)
# 
# sys.exit()

#print(outd)

denses = np.unique(d[:,1])
intensities = np.unique(d[:,2])
temps = np.unique(d[:,3])
surfs = 10**np.linspace(18.477,25.,50)

nd = denses.size
nt = temps.size
ni = intensities.size
nn = surfs.size

f = open("shrunk_table_labels_"+time.strftime("%d%m%y")+".dat",'w')

for ar in [denses,temps,intensities,surfs]:
    #outd = np.hstack([outd,ar.size,ar])
    f.write(str(ar.size)+"\n")
    for v in ar:
        f.write(str(v)+"\n")
    

f.close()

tab_slices = np.empty((nd,nt,ni),dtype=object)
# heats = np.empty((nd,nt,ni),dtype=object)
# cools = np.empty((nd,nt,ni),dtype=object)
# dusts = np.empty((nd,nt,ni),dtype=object)
# arads = np.empty((nd,nt,ni),dtype=object)
# dg = np.empty((nd,nt,ni),dtype=object)

#outp_cols_interp = [5,6,7,8,9]
outp_cols_interp = [6,7,5,8,9,10,11]
outp_dolog = [True,True,False,True,False]
noutp = len(outp_cols_interp)

alloutp = [np.empty((nd,nt,ni),dtype=object) for x in range(noutp)]

for id in range(nd):
    for it in range(nt):
        print(id,it)
        for ii in range(ni):
            tab_slices[id,it,ii] = (d[:,1]==denses[id]) & (d[:,2]==intensities[ii]) & (d[:,3]==temps[it])
            if ( np.sum(tab_slices[id,it,ii])>0 ):
                for i,icol in enumerate(outp_cols_interp):
                    alloutp[i][id,it,ii] = np.interp(surfs,d[tab_slices[id,it,ii],4],d[tab_slices[id,it,ii],icol])
#                 heats[id,it,ii] = np.interp(surfs,d[tab_slices[id,it,ii],4],d[tab_slices[id,it,ii],5])
#                 cools[id,it,ii] = np.interp(surfs,d[tab_slices[id,it,ii],4],d[tab_slices[id,it,ii],6])
#                 dusts[id,it,ii] = np.interp(surfs,d[tab_slices[id,it,ii],4],d[tab_slices[id,it,ii],7])
#                 arads[id,it,ii] = np.interp(surfs,d[tab_slices[id,it,ii],4],d[tab_slices[id,it,ii],8])
#                 dg[id,it,ii] = np.interp(surfs,d[tab_slices[id,it,ii],4],d[tab_slices[id,it,ii],9])
#                 # temporary hack for negative cooling
#                 if ( np.any(cools[id,it,ii]<0.) ):
#                     belowzero = (cools[id,it,ii]<0.)
#                     heats[id,it,ii][belowzero]-=cools[id,it,ii][belowzero]
#                     cools[id,it,ii][belowzero] = 1.e-50 # basically zero
#                 if ( np.any(heats[id,it,ii]<0.) ):
#                     belowzero = (heats[id,it,ii]<0.)
#                     cools[id,it,ii][belowzero]-=heats[id,it,ii][belowzero]
#                     heats[id,it,ii][belowzero] = 1.e-50 # basically zero


alloutdata = np.empty((noutp,0))

for id in range(nd):
    for it in range(nt):
        print(id,it)
        for ii in range(ni):
            if ( np.sum(tab_slices[id,it,ii])>0 ):
                outdata = []
                for iout,outp in enumerate(alloutp):
                    if ( outp_dolog[iout] ):
                        outdata.append(np.log10(outp[id,it,ii]))
                    else:
                        outdata.append(outp[id,it,ii])
                outdata = np.array(outdata)
                #outdata = np.array([np.log10(heats[id,it,ii]),np.log10(cools[id,it,ii]),dusts[id,it,ii],np.log10(arads[id,it,ii]),dg[id,it,ii]])
            else:
                outdata = np.array([ [0]*nn,[0]*nn,[0]*nn,[0]*nn,[0]*nn ])
            alloutdata=np.hstack([alloutdata,outdata])
                #np.savetxt(outfile,outdata,fmt='%.8e')

alloutdata = np.array(alloutdata)
np.savetxt("shrunk_table_"+time.strftime("%d%m%y")+".dat",alloutdata.T)
