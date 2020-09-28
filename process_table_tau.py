import numpy as np
import time
import sys
import pandas as pd

# loadedFile = ""
# loadedLineFile = ""
# d_in = None

class axes_offsets:
    max_error_factor = 1.e-8

    def __init__(self,d):
        self.denses = np.unique(d['n0'])
        self.intensities = np.unique(d['i0'])
        self.temps = np.unique(d['tgas'])
        self.nd = self.denses.size
        self.nt = self.temps.size
        self.ni = self.intensities.size
        self.nc = len(d)//self.nd//self.nt//self.ni
        self.d_offset = self.nt*self.ni*self.nc*np.arange(self.nd)
        self.i_offset = self.nt*self.nc*(np.arange(self.intensities.size,0,-1)-1) # intensities *decrease* through the array
        self.t_offset = self.nc*np.arange(self.temps.size)
        
        print(self.nd,self.nt,self.ni,self.nc)

    def inrange(self,id,ii,it):
        if id<0 or id>=self.nd: return False
        if ii<0 or ii>=self.ni: return False
        if it<0 or it>=self.nt: return False
        return True
    
    def convert_indices(self,other,id_in,ii_in,it_in):
        dens = other.denses[id_in]
        intensity = other.intensities[ii_in]
        temp = other.temps[it_in]
        
        for val,vals in zip([dens,intensity,temp],[self.denses,self.intensities,self.temps]):
            if val<vals.min()-self.max_error_factor or val>vals.max()+self.max_error_factor:
                return None,None,None
        
        id_out = np.abs(self.denses - dens).argmin() 
        ii_out = np.abs(self.intensities - intensity).argmin() 
        it_out = np.abs(self.temps - temp).argmin() 
        return id_out,ii_out,it_out
        

def process_table_tau(taumode,nodustmode,highdensemode,tableFile,lineTableFile=None,taumax=5.):
    if lineTableFile is not None:
            if isinstance(lineTableFile,list):
                print("Loading in line files ",lineTableFile)
                lineData = [pd.read_csv(file,delim_whitespace=True) for file in lineTableFile]
                lineMode = len(lineData)
            else:
                loadedLineFile = lineTableFile
                print("Loading in line file "+lineTableFile)
                lineData = [pd.read_csv(lineTableFile,delim_whitespace=True)]
                lineMode = 1
    else:
        lineMode = 0

    print("Loading in table file ",tableFile)
    d = pd.read_csv(tableFile,delim_whitespace=True)

    if nodustmode:
        # add in blank data
        d['tgrain'] = 0
        d['dg'] = 0

    #convert from exp(-tau) to tau
    #d[:,12] = -np.log(d[:,12]) # should already be done now # no, it's not logged, but that's fine.
    main_table_offsets = axes_offsets(d)
    
    if lineMode:
        line_table_offsets = [axes_offsets(l) for l in lineData]

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

    print("Processing")

    if taumode:
#         outp_cols_interp = [6,7,5,8,9,10,11,4]
        outp_cols_interp = ["heat","cool","tgrain","prad","dg","kabs","kscat","colden"]
        icol_in = "tau"
        outp_dolog = [True,True,False,True,False,False,False,False]
    else:
#         outp_cols_interp = [6,7,5,8,9,10,11,12]
        outp_cols_interp = ["heat","cool","tgrain","prad","dg","kabs","kscat","tau"]
        icol_in = "colden"
        outp_dolog = [True,True,False,True,False,False,False,False]
    noutp = len(outp_cols_interp)
    
#     nlines=7
    linecols = ["CO1","CO2","HCN1","HCN2","H2_1","H2_2","H2_3","nujnu12","nujnu8","nujnu850"]
    nlines=len(linecols)
    noutp+=nlines # for lines
    outp_dolog+=[False]*nlines
    line_coldens_col = "colden"

    alloutp = [np.empty((main_table_offsets.nd,main_table_offsets.nt,main_table_offsets.ni),dtype=object) for x in range(noutp)]
    

#     nc = d.shape[0]//nd//nt//ni
#     d_offset = nt*ni*nc*np.arange(nd)
#     i_offset = nt*nc*(np.arange(intensities.size,0,-1)-1) # intensities *decrease* through the array
#     t_offset = nc*np.arange(temps.size)

    print("Interpolating")
    for id in range(main_table_offsets.nd):
        for it in range(main_table_offsets.nt):
#             print(id,it)
            for ii in range(main_table_offsets.ni):
                index0 = main_table_offsets.d_offset[id]+main_table_offsets.i_offset[ii]+main_table_offsets.t_offset[it]
                index1 = index0+main_table_offsets.nc

                for i,icol in enumerate(outp_cols_interp):
                    alloutp[i][id,it,ii] = np.interp(column_in,d[icol_in][index0:index1],d[icol][index0:index1])
#                     if icol in d:
#                         alloutp[i][id,it,ii] = np.interp(column_in,d[icol_in][index0:index1],d[icol][index0:index1])
#                     else:
#                         alloutp[i][id,it,ii] = np.zeros((ncol_in))


                for i in range(noutp-nlines,noutp):
                    alloutp[i][id,it,ii] = np.zeros((ncol_in)) # initalize to 0

                if lineMode:
                    for itable,(ld,lto) in enumerate(zip(lineData,line_table_offsets)):
                        l_id,l_ii,l_it = lto.convert_indices(main_table_offsets,id,ii,it)
                        if l_id is not None:
                            index0 = lto.d_offset[l_id]+lto.i_offset[l_ii]+main_table_offsets.t_offset[l_it]
                            index1 = index0+lto.nc
                            if taumode:
                                col_select = np.interp(column_in,d[icol_in][index0:index1],d['colden'][index0:index1])
        #                         print(taumode,col_select.shape)
                                for i,icol in enumerate(linecols): 
                                    if icol in ld:
                                        alloutp[i+noutp-nlines][id,it,ii] = np.interp(col_select,ld[line_coldens_col][index0:index1],ld[icol][index0:index1])
                            else:
        #                         print(taumode,column_in.shape)
                                for i,icol in enumerate(linecols):
                                    if icol in ld:
                                        alloutp[i+noutp-nlines][id,it,ii] = np.interp(column_in,ld[line_coldens_col][index0:index1],ld[icol][index0:index1])

    print("Collating")
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
