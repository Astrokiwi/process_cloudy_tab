import process_table_tau
import h_weighted_table
import time
from multiprocessing import Pool
from functools import partial
import sys

import numpy as np

# masses = [1.e-11,1.e-10,1.e-9,1.e-8,1.e-7,1.e-6,1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.,10.]
# masses = [1.e-4]


# masses = [0.1] # for large scale production runs
# masses = [1.e-4] # for high res production runs

masses = [1.e-4,1.e-3,1.e-2,1.e-1,1.,10.,100.] # res scale run for large scale production runs

# masses = np.linspace(0.01,1.,20)
# masses = [0.01,0.05,0.1]

# masses = [1.e-4]

stochastic_accel = True

# full for production
table_parameters = [ [False,True,"highden_260118.txt"],
                     [False,False,"tables_100818.txt",["linetables_281118.txt","emissivity_170620.txt"]],
                     [True,False,"nodust_301117.txt"]
                                            ]

# quick for test
# table_parameters = [ 
#                      [False,False,"tables_100818.txt",["linetables_281118.txt","emissivity_170620.txt"]],
#                                             ]

# table_parameters = [ [False,False,"tables_100818.txt","linetables_281118.txt"]
#                                             ]

# table_parameters = [ 
#                      [False,False,"tables_100818.txt","linetables_091118.txt"]
#                                             ]

table_parameters_taumode = [ [True]+line for line in table_parameters ] + [ [False]+line for line in table_parameters ]



table_parameters_mass = []
for mass in masses:
    table_parameters_mass = table_parameters_mass + [line[0:3]+[mass] for line in table_parameters]

table_date = time.strftime("%d%m%y")
for line in table_parameters_mass:
    line[2] = table_date

def process_table_tab_for_pool(prams):
    process_table_tau.process_table_tau(*prams,taumax=7.)

def generate_h_weighted_table_for_pool(prams):
    h_weighted_table.generate_h_weighted_table(*prams,stochastic_accel=stochastic_accel)

with Pool(processes=64) as pool:
    pool.map(process_table_tab_for_pool,table_parameters_taumode)
    pool.map(generate_h_weighted_table_for_pool,table_parameters_mass)

