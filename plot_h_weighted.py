import numpy as np

from sys import path,exit
path.append("../src/")
import tab_interp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

import itertools
marker = itertools.cycle((',', '+', '.', 'o', '*')) 
color_cycle = [x[u'color'] for x in mpl.rcParams['axes.prop_cycle']]

print("Loading tables")

# mass_p =[1.e-5,5.e-5,1.e-4,5.e-4,1.e-3,5.e-3,1.e-2]
# mass_p =[1.e-5,1.e-4,1.e-3,1.e-2,3.e-2,1.e-1,3.e-1,1.,10.]
# mass_p =[1.e-11,1.e-10,1.e-9,1.e-8,1.e-7,1.e-6,1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.,10.,100.]
mass_p =[1.e-4,1.e-3,1.e-2,1.e-1,1.,10.,100.]

# mass_p = np.linspace(0.01,1.,20)
# mass_p = [0.01,0.05,0.1]

# mass_p =[1.e-5,5.e-5,1.e-4,5.e-4,1.e-3,5.e-3,1.e-2]


# tau_tables = [tab_interp.CoolHeatTab("shrunk_table_labels_171117tau.dat","shrunk_table_281117_m{}_hsmooth_tau.dat".format(x)) for x in mass_p]
# tau_tables = [tab_interp.CoolHeatTab( ("shrunk_table_labels_291117tau.dat"),
#                                             ("shrunk_table_291117_m{}_hsmooth_tau.dat".format(x)),
#                                             ("shrunk_table_labels_011217taunodust.dat"),
#                                             ("shrunk_table_011217_m0.1_hsmooth_taunodust.dat")
#                                             ) for x in mass_p]

# date = "181018"
date = "010519"

tau_tables = [tab_interp.CoolHeatTab( ("shrunk_table_labels_{}tau.dat".format(date)),
                                        ("shrunk_table_{}_m{}_hsmooth_tau.dat".format(date,x)),
                                        ("shrunk_table_labels_{}taunodust.dat".format(date)),
                                        ("shrunk_table_{}_m{}_hsmooth_taunodust.dat".format(date,x)),
                                        ("shrunk_table_labels_{}taudense.dat".format(date)),
                                        ("shrunk_table_{}_m{}_hsmooth_taudense.dat".format(date,x))
                                        ) for x in mass_p]

table_attributes = ['dHeat','dCool','dustT','arad']#,'dg','opac_abs','opac_scat','column_out','line_co1','line_co2','line_hcn1','line_hcn2']

# P.rc('text', usetex=True)


for tau_p in [0.25,0.5,0.75]:
    for logT in [2,3,4]:
# for tau_p in [0]:
#     for logT in [2]:

        fig,ax1=P.subplots()
        ax1.set_xlabel(r"$\log M_p$ (M$_\odot$)")
        # ax1.set_ylabel(r"$\log a_r$ (cm s$^{-2}$)")
        ax1.set_ylabel(r"$a_r$ (cm s$^{-2}$)")
        maxy = 0.
        for lognH in [1,2,3,4,5,6]:
            nH_p = 10.**lognH

        #     nH_p = 10**6.
            TK_p = 10.**logT
            flux_p = 10.**6.
#             tau_p = 0.

            def m_to_h(m):
                return nH_p**(-1./3.)/0.3404*(m/0.1)**(1./3.)

            print("Interpolating from tables")

            interp_structs = [table.interpTab(nH_p,TK_p,flux_p,tau_p) for table in tau_tables]

            vals = [[struct.__getattr__(attribute) for attribute in table_attributes] for struct in interp_structs]

            vals = np.array(vals)

            vals = np.vstack((vals.T,mass_p)).T

            # np.savetxt("../data/h_weighted_test.txt",vals)

        #     ax1.scatter(np.log10(mass_p),vals[:,3])
            a_p = 10.**vals[:,3]
            
            maxy = max(np.max(a_p),maxy)
#             print(vals[:,3],a_p)
            lastmarker = next(marker)
            ax1.scatter(np.log10(mass_p),a_p,label=r"$n_H=10^%d$ cm$^{-3}$"%lognH,marker=lastmarker)
            if lognH==4:
                ax2 = ax1.twiny()
                #     ax2.scatter(np.log10(mass_p),vals[:,3])
                n = len(mass_p)
                transparent_colours = np.array([(1.,1.,1.,0.)]*n)
                ax2.scatter(np.log10(mass_p),a_p,c=transparent_colours)
                masslocs = ax2.get_xticks()
                m = 10.**masslocs
                h_p = m_to_h(m)
#                 print(m)
#                 print(h_p)
                ax2.set_xticklabels(["%.1f"%x for x in np.log10(h_p)])
                ax2.set_xlabel(r"$\log h_p$ (pc) ($n_H=10^%d$ cm$^{-3}$)"%lognH)

        ax1.set_ylim([0,maxy])
        ax1.legend()

            # P.yscale('log')
            # P.xscale('log')



        P.savefig("../../figures/h_weighted_test_T{}_tau{}.pdf".format(logT,tau_p))
        P.close('all')
