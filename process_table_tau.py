import numpy as np
import time
import sys

# loadedFile = ""
# loadedLineFile = ""
# d_in = None

class axes_offsets:
    def __init__(self,d):
        self.denses = np.unique(d[:,1])
        self.intensities = np.unique(d[:,2])
        self.temps = np.unique(d[:,3])
        self.nd = self.denses.size
        self.nt = self.temps.size
        self.ni = self.intensities.size
        self.nc = d.shape[0]//self.nd//self.nt//self.ni
        self.d_offset = self.nt*self.ni*self.nc*np.arange(self.nd)
        self.i_offset = self.nt*self.nc*(np.arange(self.intensities.size,0,-1)-1) # intensities *decrease* through the array
        self.t_offset = self.nc*np.arange(self.temps.size)
        
        print(self.nd,self.nt,self.ni,self.nc)

    def inrange(self,id,ii,it):
        if id<0 or id>=self.nd: return False
        if ii<0 or ii>=self.ni: return False
        if it<0 or it>=self.nt: return False
        return True

def process_table_tau(taumode,nodustmode,highdensemode,tableFile,lineTableFile=None,taumax=5.):
    global loadedFile,d_in,loadedLineFile,lineData
    
    if 'loadedFile' in globals() and tableFile==loadedFile:
        print("Using table already in memory - hopefully you want to do this!")
    else:
        print("Loading in table file ",tableFile)
        loadedFile = tableFile
        d_in = np.loadtxt(tableFile,skiprows=1)
    
    if lineTableFile is not None:
        if 'loadedLineFile' in globals() and loadedLineFile==lineTableFile:
            print("Using line table already in memory - hopefully you want to do this!")
            lineMode = True
        else:
            loadedLineFile = lineTableFile
            print("Loading in line file "+lineTableFile)
            lineData = np.loadtxt(lineTableFile,skiprows=1)
            lineMode = True
    else:
        lineMode = False

# if re-running in same name-space, don't need to reload the data
# if ( not 'd_in' in dir() ):
#     print("Loading in table")
#     # load in giant file
#     d_in = np.loadtxt("highden_260118.txt",skiprows=1)
#     #d_in = np.loadtxt("tables_271117.txt",skiprows=1)
# #     d_in = np.loadtxt("nodust_301117.txt",skiprows=1)
# #     d_in = np.loadtxt("highden_260118.txt",skiprows=1)
# #     d_in = np.loadtxt("tables_271117.txt",skiprows=1)
# #     d_in = np.loadtxt("nodust_301117.txt",skiprows=1)
# else:
#     print("Using table already in memory - hopefully you want to do this!")

    d = np.copy(d_in)

    if nodustmode:
        # add in blank data
        s = d.shape[0]
        d = np.insert(d,5,np.zeros(s),axis=1)
        d = np.insert(d,9,np.zeros(s),axis=1)

    #convert from exp(-tau) to tau
    #d[:,12] = -np.log(d[:,12]) # should already be done now # no, it's not logged, but that's fine.
    main_table_offsets = axes_offsets(d)
    
    if lineMode:
        line_table_offsets = axes_offsets(lineData)

#     denses = np.unique(d[:,1])
#     intensities = np.unique(d[:,2])
#     temps = np.unique(d[:,3])

    if taumode:
        column_in = 10.**np.linspace(-2.,np.log10(taumax),50)
    else:
        column_in = 10.**np.linspace(18.477,25.,50)
    ncol_in = column_in.size

    #taus = np.linspace(0.,5.,50)

#     nd = denses.size
#     nt = temps.size
#     ni = intensities.size

    tau_suffix = "tau" if taumode else "coldens"
    dust_suffix = "nodust" if nodustmode else ""
    if nodustmode and highdensemode:
        raise Exception("Can't have both nodustmode and highdensemode!")
    dense_suffix = "dense" if highdensemode else ""
 
    suffixes = tau_suffix+dust_suffix+dense_suffix
    print("Saving labels")
    f = open("shrunk_table_labels_"+time.strftime("%d%m%y")+suffixes+".dat",'w')

    for ar in [main_table_offsets.denses,main_table_offsets.temps,main_table_offsets.intensities,column_in]:
        #outd = np.hstack([outd,ar.size,ar])
        f.write(str(ar.size)+"\n")
        for v in ar:
            f.write(str(v)+"\n")
    

    f.close()

    if taumode:
        outp_cols_interp = [6,7,5,8,9,10,11,4]
        icol_in = 12
        outp_dolog = [True,True,False,True,False,False,False,False,False,False,False,False]
    else:
        outp_cols_interp = [6,7,5,8,9,10,11,12]
        icol_in = 4
        outp_dolog = [True,True,False,True,False,False,False,False,False,False,False,False]
    noutp = len(outp_cols_interp)
    
    noutp+=4 # for lines
    linecols = [5,6,7,8]
    line_coldens_col = 4

    alloutp = [np.empty((main_table_offsets.nd,main_table_offsets.nt,main_table_offsets.ni),dtype=object) for x in range(noutp)]

#     nc = d.shape[0]//nd//nt//ni
#     d_offset = nt*ni*nc*np.arange(nd)
#     i_offset = nt*nc*(np.arange(intensities.size,0,-1)-1) # intensities *decrease* through the array
#     t_offset = nc*np.arange(temps.size)


    for id in range(main_table_offsets.nd):
        for it in range(main_table_offsets.nt):
#             print(id,it)
            for ii in range(main_table_offsets.ni):
                index0 = main_table_offsets.d_offset[id]+main_table_offsets.i_offset[ii]+main_table_offsets.t_offset[it]
                index1 = index0+main_table_offsets.nc
#                 d_slice = d[index0:index1]

                for i,icol in enumerate(outp_cols_interp):
                    alloutp[i][id,it,ii] = np.interp(column_in,d[index0:index1,icol_in],d[index0:index1,icol])
                if lineMode and line_table_offsets.inrange(id,ii,it):
                    index0 = line_table_offsets.d_offset[id]+line_table_offsets.i_offset[ii]+main_table_offsets.t_offset[it]
                    index1 = index0+line_table_offsets.nc
                    if taumode:
                        col_select = np.interp(column_in,d[index0:index1,icol_in],d[index0:index1,4])
#                         print(taumode,col_select.shape)
                        for i,icol in enumerate(linecols):
                            alloutp[i+noutp-4][id,it,ii] = np.interp(col_select,lineData[index0:index1,line_coldens_col],lineData[index0:index1,icol])
                    else:
#                         print(taumode,column_in.shape)
                        for i,icol in enumerate(linecols):
                            alloutp[i+noutp-4][id,it,ii] = np.interp(column_in,lineData[index0:index1,line_coldens_col],lineData[index0:index1,icol])
                else:
#                     print(taumode,lineMode,ncol_in)
                    for i in range(noutp-4,noutp):
                        alloutp[i][id,it,ii] = np.zeros((ncol_in)) # placeholder

    alloutdata = np.empty((noutp,0))

    for id in range(main_table_offsets.nd):
        for it in range(main_table_offsets.nt):
#             print(id,it)
            for ii in range(main_table_offsets.ni):
                outdata = []
                for iout,outp in enumerate(alloutp):
                    if ( outp_dolog[iout] ):
                        outdata.append(np.log10(outp[id,it,ii]))
                    else:
                        outdata.append(outp[id,it,ii])
                outdata = np.array(outdata)
#                 print(outdata.shape,alloutdata.shape)
                alloutdata=np.hstack([alloutdata,outdata])

    alloutdata = np.array(alloutdata)
    print("Saving table")
    np.savetxt("shrunk_table_"+time.strftime("%d%m%y")+suffixes+".dat",alloutdata.T)

if __name__ == '__main__':

    #taumode=True means use tau as the independent variable
    #taumode=False means use NH column density instead
    taumode = False

    #nodustmode = this is dust-free gas for the secondary table
    nodustmode = False

    #highdensemode = this is the third table, high density but only cold gas
    highdensemode = False

    process_table_tau(taumode,nodustmode,highdensemode,"tables_100818.txt","linetables_091118.txt")
