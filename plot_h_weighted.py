import numpy as np

from sys import path,exit
path.append("../src/")
import tab_interp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P

print("Loading tables")

# mass_p =[1.e-5,5.e-5,1.e-4,5.e-4,1.e-3,5.e-3,1.e-2]
# mass_p =[1.e-5,1.e-4,1.e-3,1.e-2,3.e-2,1.e-1,3.e-1,1.,10.]
# mass_p =[1.e-5,1.e-4,1.e-3,1.e-2,3.e-2,1.e-1,3.e-1,1.,10.]
# mass_p = np.linspace(0.01,1.,20)
# mass_p = [0.01,0.05,0.1]

mass_p =[1.e-5,5.e-5,1.e-4,5.e-4,1.e-3,5.e-3,1.e-2]


# tau_tables = [tab_interp.CoolHeatTab("shrunk_table_labels_171117tau.dat","shrunk_table_281117_m{}_hsmooth_tau.dat".format(x)) for x in mass_p]
# tau_tables = [tab_interp.CoolHeatTab( ("shrunk_table_labels_291117tau.dat"),
#                                             ("shrunk_table_291117_m{}_hsmooth_tau.dat".format(x)),
#                                             ("shrunk_table_labels_011217taunodust.dat"),
#                                             ("shrunk_table_011217_m0.1_hsmooth_taunodust.dat")
#                                             ) for x in mass_p]

date = "140818"

tau_tables = [tab_interp.CoolHeatTab( ("shrunk_table_labels_{}tau.dat".format(date)),
                                        ("shrunk_table_{}_m{}_hsmooth_tau.dat".format(date,x)),
                                        ("shrunk_table_labels_{}taunodust.dat".format(date)),
                                        ("shrunk_table_{}_m{}_hsmooth_taunodust.dat".format(date,x)),
                                        ("shrunk_table_labels_{}taudense.dat".format(date)),
                                        ("shrunk_table_{}_m{}_hsmooth_taudense.dat".format(date,x))
                                        ) for x in mass_p]

table_attributes = ['dHeat','dCool','dustT','arad','dg','opac_abs','opac_scat','column_out']

nH_p = 10**6.
TK_p = 10.**2.5
flux_p = 10.**6.
tau_p = 0.

def m_to_h(m):
    return nH_p**(-1./3.)/0.3404*(m/0.1)**(1./3.)

print("Interpolating from tables")

interp_structs = [table.interpTab(nH_p,TK_p,flux_p,tau_p) for table in tau_tables]

vals = [[struct.__getattr__(attribute) for attribute in table_attributes] for struct in interp_structs]

vals = np.array(vals)

vals = np.vstack((vals.T,mass_p)).T

# np.savetxt("../data/h_weighted_test.txt",vals)

fig,ax1=P.subplots()
ax1.scatter(np.log10(mass_p),vals[:,3])
# ax1.scatter(mass_p,vals[:,3])
ax1.set_xlabel(r"$\log M_p$ (M$_\odot$)")
ax1.set_ylabel(r"$\log a_r$ (cm s$^{-2}$)")

# ax2 = ax1.twiny()
# ax2.scatter(np.log10(mass_p),vals[:,3])
# # P.tight_layout()
# masslocs = ax2.get_xticks()
# m = 10.**masslocs
# h_p = m_to_h(m)
# ax2.set_xticklabels(["%.1f"%x for x in np.log10(h_p)])
# ax2.set_xlabel(r"$\log h_p$ (pc)")

# P.yscale('log')
# P.xscale('log')
P.savefig("../../figures/h_weighted_test.pdf")
P.close('all')
