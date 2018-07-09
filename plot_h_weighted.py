import numpy as np

from sys import path,exit
path.append("../src/")
import tab_interp

import matplotlib.pyplot as P

print("Loading tables")

mass_p = np.arange(1,11)*0.01
#mass_p = np.array([0.01,0.02])

# mass_p = [0.02,0.04,0.08,0.1]

# tau_tables = [tab_interp.CoolHeatTab("shrunk_table_labels_171117tau.dat","shrunk_table_281117_m{}_hsmooth_tau.dat".format(x)) for x in mass_p]
tau_tables = [tab_interp.CoolHeatTab( ("shrunk_table_labels_291117tau.dat"),
                                            ("shrunk_table_291117_m{}_hsmooth_tau.dat".format(x)),
                                            ("shrunk_table_labels_011217taunodust.dat"),
                                            ("shrunk_table_011217_m0.1_hsmooth_taunodust.dat")
                                            ) for x in mass_p]

table_attributes = ['dHeat','dCool','dustT','arad','dg','opac_abs','opac_scat','column_out']

nH_p = 10**6.
TK_p = 10.**2.5
flux_p = 10.**6.
tau_p = 0.

print("Interpolating from tables")

interp_structs = [table.interpTab(nH_p,TK_p,flux_p,tau_p) for table in tau_tables]

vals = [[struct.__getattr__(attribute) for attribute in table_attributes] for struct in interp_structs]

vals = np.array(vals)

vals = np.vstack((vals.T,mass_p)).T

#np.savetxt("h_weighted_test.txt",vals)

P.figure()
P.scatter(mass_p,vals[:,3])
P.xlabel(r"$M_p$ (M$_\odot$)")
P.ylabel(r"$\log a_r$ (cm s$^{-2}$)")
P.savefig("../../figures/h_weighted_test.pdf")
P.close('all')
