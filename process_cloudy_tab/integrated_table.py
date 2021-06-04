# by integrating 1/(heating-cooling), we can do more exact heating/cooling calculations without needing a short timestep
# requires python 2 for some reason?

import numpy as np
import time

import pandas
import scipy.integrate as integrate

table_date = "140618"
table_labels_date = "240518"

table_axis_labels = ['dense','temp','intensity','tau']
naxes = 4

molecular_mass = 4./(1.+3.*.76)
proton_mass_cgs = 1.6726e-24
gamma_minus_one = 5./3.-1.
boltzmann_cgs = 1.38066e-16

U_TO_TK = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)

for dust_suffix in ["","nodust"]:
# for dust_suffix in [""]:
    labelFileName = "shrunk_table_labels_{}tau{}.dat".format(table_labels_date,dust_suffix)
    tableFileName = "shrunk_table_{}_m0.01_hsmooth_tau{}.dat".format(table_date,dust_suffix)
    labelVals = np.loadtxt(labelFileName)
    
    offset = 0
    axis_l = dict()
    axis_points = dict()
    for axis_label in table_axis_labels:
        axis_l[axis_label] = int(labelVals[offset])
        offset+=1
        axis_points[axis_label] = labelVals[offset:offset+axis_l[axis_label]]
        offset+=axis_l[axis_label]
    
    for axis_label in table_axis_labels:
        print(axis_label,axis_l[axis_label],axis_points[axis_label])

    linearTempPoints = 10.**axis_points['temp']
    massDensePoints = 10.**axis_points['dense']*molecular_mass*proton_mass_cgs
    print(massDensePoints)
    denseStepSize = axis_l['tau']*axis_l['temp']*axis_l['intensity']
    
    
    tableData = pandas.read_csv(tableFileName,delim_whitespace=True,header=None,names=['dHeat','dCool','dustT','arad','dg','opac_abs','opac_scat','coldens'])
    allDense = np.repeat(massDensePoints,denseStepSize)
    tableData['inverseHeatCool'] = allDense/(10.**tableData['dHeat']-10.**tableData['dCool'])
    tableData['integratedInverseHeatCool'] = tableData['dHeat']*0.
    tableData['opac_sum'] = tableData['opac_abs']+tableData['opac_scat']
    print(tableData.columns.values)
    tableData=tableData.drop(['dHeat','dCool','opac_abs','opac_scat','coldens'],1)
    print(tableData.columns.values)
    
    print(np.max(tableData['inverseHeatCool']))
#     tableData['heatCool'] = (10.**tableData['dHeat']-10.**tableData['dCool'])
    
    linearUPoints = linearTempPoints/U_TO_TK
    
    # test - just one set of densities etc
    
    for i_dense in range(axis_l['dense']):
        print(i_dense)
        for i_intense in range(axis_l['intensity']):
            for i_tau in range(axis_l['tau']):
                offset = i_tau + i_intense*axis_l['tau'] + i_dense*axis_l['tau']*axis_l['temp']*axis_l['intensity']
                istep = axis_l['tau']*axis_l['intensity']
                heating_slice = tableData['inverseHeatCool'][offset:offset+istep*axis_l['temp']:istep]
            #     print(np.array(zip(axis_points['temp'],heating_slice)))
                integ = integrate.cumtrapz(heating_slice,x=linearUPoints,initial=0.)
                tableData['integratedInverseHeatCool'][offset:offset+istep*axis_l['temp']:istep] = integ
            #     print(np.array(zip(axis_points['temp'],integ-integ[2])))
            #     print(tableData['inverseHeatCool'])


    tableData=tableData.drop('inverseHeatCool',1)
    print(np.min(tableData['integratedInverseHeatCool']),np.max(tableData['integratedInverseHeatCool']))
    
    # xxxx don't do this, it makes the fit horrible later
    # add offset so that integrated heatcool is never below zero
#     minheatcool = np.min(tableData['integratedInverseHeatCool'])
#     if minheatcool<0.:
#         tableData['integratedInverseHeatCool']-=minheatcool
    

    tableData.to_csv("integ_table_{}_m0.01_hsmooth_tau{}.dat".format(table_date,dust_suffix),sep=' ')






















