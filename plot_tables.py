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
    d_in = np.loadtxt("tables_251017.txt",skiprows=1)
else:
    print("Using table already in memory - hopefully you want to do this!")

class Data_Axis:
    def __init__(self,n,vals,i_dti,str):
        self.n = n
        self.vals = vals
        self.i_dti = i_dti
        self.str = str
    
    def dti(self,id,it,ii):
        l = [id,it,ii]
        return l[self.i_dti]
    
    def label_str(self,i):
        return self.str+"{}".format(self.vals[i])
    
print("Working on table")

#cutoffs = [23,24,25,26]
#cutoffs = [24]
cutoffs = [2.]

for cutoff in cutoffs:

    d = np.copy(d_in) # so we don't repeat the logging, but we can change the logging without having to read in again

    # use tau for x axis
    d[:,4] = -np.log(d[:,12])

    # set bound (optional)
    #coldens_cutoff = 1.e25
    #coldens_cutoff = 10.**cutoff
#     coldens_cutoff = cutoff
#     coldens_slice = (d[:,4]<coldens_cutoff)
#     d = d[coldens_slice,:]

    # coldens_cutoff = 1.e16
    # coldens_slice = (d[:,4]>coldens_cutoff)
    # d = d[coldens_slice,:]

    # create new data (optional)
    # x = d[:,12] * (d[:,10]+d[:,11]) * 10.**d[:,2] # approximate rad pressure modulo units = exp(-tau) * (sum of opacities) * input intensity 
    # d = np.vstack((d.T,x)).T # add to array - 14th column
    # d[:,13] = np.log10(d[:,13]) # and log this

    d[:,12] = -np.log(d[:,12]) # convert from exp(-tau) to tau
    #d[:,12] = -np.log(np.clip(d[:,12],None,1.+1.e-10)) # convert from exp(-tau) to tau
    # d[:,12]/=d[:,4] # divide by column density to get more "visibility" at low column density - doesn't work, just gives flat stuff :/

    #d[:,4:9] = np.log10(d[:,4:9]) # log most dependent variables, plus column density (other independent variables are already logged) (don't log dust fraction)
    d[:,5:9] = np.log10(d[:,5:9]) # log most dependent variables, BUT NOT column density (other independent variables are already logged) (don't log dust fraction)
    d[:,10:12] = np.log10(d[:,10:12]) # log opacities too
    #d[:,12] = np.log10(d[:,12]) # log tau too


    d[:,6]-=d[:,1] # convert from erg/s/cm^3 to erg/s/(particle)
    d[:,7]-=d[:,1] # convert from erg/s/cm^3 to erg/s/(particle)

    # find all values for independent variables (except column density)
    denses = np.unique(d[:,1])
    intensities = np.unique(d[:,2])
    temps = np.unique(d[:,3])

    nd = denses.size
    nt = temps.size
    ni = intensities.size

    nc = d.shape[0]/nd/nt/ni
    
#     tau0 = d[:nc,12].copy()
#     for i in range(nt*ni*nd):
#         d[i*nc:(i+1)*nc,12]/=tau0
#     d[:,12] = np.log10(d[:,12])
#     bads = ~np.isfinite(d[:,12])
#     d[bads,12] = np.nan


    data_axes = dict()
    data_axes["d"] = Data_Axis(nd,denses,0,"n")
    data_axes["t"] = Data_Axis(nt,temps,1,"T")
    data_axes["i"] = Data_Axis(ni,intensities,2,"I")

    # output all different quantities
    titles = ["tgrain","heat","cool","prad","dg","kabs","kscat","tau"]
    #toplot = [False,False,False,True,False,False,False,False]
    #toplot = [False,False,False,False,False,False,False,True]
    toplot = [True]*8
    #titles = ["tgrain","heat"]
    #titles = ["tgrain","heat","cool","prad","dg","kabs","kscat","tau","effec_prad"]

    #colrowline_var0 = np.array(("t","d","i"))
    #colrowline_var_list = [colrowline_var0,np.roll(colrowline_var0,1),np.roll(np.roll(colrowline_var0,1),1)]
    #colrowline_var_list = [("t","d","i")]
    colrowline_var_list = [("t","i","d")]

    for colrowline_var in colrowline_var_list:

    ## vary across columns of plots
    #column_var = "i"
    ## vary across rows of plots
    #row_var = "d"
    ## vary within each plot
    #line_var = "t"

        # vary across columns of plots
        column_var = colrowline_var[0]
        # vary across rows of plots
        row_var = colrowline_var[1]
        # vary within each plot
        line_var = colrowline_var[2]

        print("Calculating order of arrays")
        nc = d.shape[0]/nd/nt/ni
        d_offset = nt*ni*nc*np.arange(nd)
        i_offset = nt*nc*(np.arange(intensities.size,0,-1)-1) # intensities *decrease* through the array
        t_offset = nc*np.arange(temps.size)


        nx = data_axes[column_var].n
        ny = data_axes[row_var].n

        figs = []
        sps = []
        yranges = [[np.nanmin(d[:,icol+5]),np.nanmax(d[:,icol+5])] for icol in range(len(titles))]
        #xrange = [np.nanmin(d[:,4]),np.nanmax(d[:,4])]
        xrange = [np.nanmin(d[:,4]),cutoff]
        x_size = xrange[1]-xrange[0]
        x_offset = x_size*.1
        y_sizes = [yr[1]-yr[0] for yr in yranges]
        y_offsets = [ys*.1 for ys in y_sizes]

        print(nx,"x",ny,"plots",colrowline_var)



        print("Setting up subplots")
        # set up figures
        for ifig, title in enumerate(titles):
            print(title)
            fig,sp = P.subplots(1,1,figsize=(2.4*nx,1.5*ny),dpi=100)
            fig.suptitle(title)
            figs.append(fig)
            sps.append(sp)

        color_cycle = [x[u'color'] for x in mpl.rcParams['axes.prop_cycle']]

        def doplot(d,sps,denses,intensities,temps,id,ni,it):
            print(id,it)
            for ii in range(ni):
                if data_axes[row_var].dti(id,it,ii)==0 and data_axes[column_var].dti(id,it,ii)==0:
                    setuplegend = True
                    leg_ax = data_axes[line_var]
                    legendlabels = map(leg_ax.label_str,range(leg_ax.n))
                    #print("making legend",legendlabels)
                else:
                    setuplegend = False
                #x_left = it*(x_size+x_offset)
                x_left = data_axes[column_var].dti(id,it,ii)*(x_size+x_offset)
            #for ii in [5]: # just one intensity
        #         tab_slice = (d[:,1]==denses[id]) & (d[:,2]==intensities[ii]) & (d[:,3]==temps[it])
                index0 = d_offset[id]+i_offset[ii]+t_offset[it]
                #index1 = index0+nc
                index1 = index0+np.searchsorted(d[index0:index0+nc,4],cutoff) # for plotting tau cutoff
                d_slice = d[index0:index1]
                #d_slice = d_slice[::10,:] # for speed
        
        #         for icol,fig in enumerate(figs):
                for icol,sp in enumerate(sps):
                    if ( toplot[icol] ):
                        #if ( icol==7 ):
                        #    d_slice[:,5+icol]/=tau0
                    
                        #y_bottom = id*(y_sizes[icol]+y_offsets[icol])
                        y_bottom = data_axes[row_var].dti(id,it,ii)*(y_sizes[icol]+y_offsets[icol])
                        color_index = data_axes[line_var].dti(id,it,ii)%len(color_cycle)
                        #color_string = 'C%1d'%color_index
                        if setuplegend: leglabel = legendlabels[leg_ax.dti(id,it,ii)]
                        else: leglabel = None
                        sp.plot(d_slice[:,4]-xrange[0]+x_left,d_slice[:,5+icol]-yranges[icol][0]+y_bottom,c=color_cycle[color_index],label=leglabel)

        print("Plotting")

        #plot all the subplots
        for id in range(nd):
            for it in range(nt):
                doplot(d,sps,denses,intensities,temps,id,ni,it)

        for ix in range(data_axes[column_var].n):
            for iy in range(data_axes[row_var].n):
                for icol,sp in enumerate(sps):
                    caption = data_axes[column_var].label_str(ix)+data_axes[row_var].label_str(iy)
                    sp.text(ix*(x_size+x_offset),(iy+1)*y_sizes[icol]+iy*y_offsets[icol],caption)
                    #sp.text(it*(x_size+x_offset),(id+1)*y_sizes[icol]+id*y_offsets[icol],"n{}T{}".format(denses[id],temps[it]))

        xticks_left = np.arange(nx)*(x_size+x_offset)
        xticks_right = np.arange(nx)*(x_size+x_offset)+x_size

        nxticks = 7
        nyticks = 5

        xticks_each = np.linspace(0,x_size,nxticks)
        xticks_all = np.tile(xticks_each,nx)+np.repeat(xticks_left,nxticks)

        #xtickvalues = d[np.linspace(0,nc-1,nxticks,dtype=int),4]
        xtickvalues = np.linspace(xrange[0],xrange[1],nxticks)
        #xticklabels = np.tile(map(lambda x : "%.1f" % x,xtickvalues),nt)
        xticklabels = np.tile(map(lambda x : "%.2e" % x,xtickvalues),nx)

        print("Rendering and dumping")
        # format & dump everything
        for ifig, fig in enumerate(figs):
            if ( toplot[ifig] ):

                sp = sps[ifig]
                sp.set_xlim([0,nx*x_size+nx*x_offset])
                sp.set_ylim([0,ny*y_sizes[ifig]+ny*y_offsets[ifig]])
            #     sp.set_xlim([0,nt*x_size+nt*x_offset])
            #     sp.set_ylim([0,nd*y_sizes[ifig]+nd*y_offsets[ifig]])
                sp.set_xticks(xticks_all)
                sp.set_xticklabels(xticklabels,rotation=90)


                yticks_bottom = np.arange(ny)*(y_sizes[ifig]+y_offsets[ifig])
                yticks_each = np.linspace(0,y_sizes[ifig],nyticks)
                yticks_all = np.tile(yticks_each,ny)+np.repeat(yticks_bottom,nyticks)

                ytickvalues = np.linspace(yranges[ifig][0],yranges[ifig][1],nyticks)
                yticklabels = np.tile(map(lambda x : "%.2e"%x,ytickvalues),ny)
                #yticklabels = np.tile(map(lambda x : "{}".format(x),ytickvalues),nd)
    
                sp.set_yticks(yticks_all)
                sp.set_yticklabels(yticklabels)

                #sp.legend()

                #handles, labels = sp.get_legend_handles_labels()
                #fig.legend(handles, labels, loc='upper left', ncol=3)

                fig.tight_layout()
                #fig.savefig("../../figures/table_summary_final_"+"".join(colrowline_var)+titles[ifig]+"{}.png".format(cutoff))
                fig.savefig("../../figures/table_summary_tau_"+"".join(colrowline_var)+titles[ifig]+"{}.png".format(cutoff))

        # to avoid memory getting full of figures
        P.close('All')
