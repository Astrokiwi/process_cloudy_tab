import process_table_tau
import h_weighted_table
import time
from multiprocessing import Pool
from functools import partial
import sys

import numpy as np

masses = [1.e-5,1.e-4,1.e-3,1.e-2,3.e-2,1.e-1,3.e-1,1.,10.]

# masses = np.linspace(0.01,1.,20)
# masses = [0.01,0.05,0.1]

# masses = [1.e-4]

table_parameters = [ [False,True,"highden_260118.txt"],
                     [False,False,"tables_100818.txt"],
                     [True,False,"nodust_301117.txt"]
                                            ]

table_parameters_taumode = [ [True]+line for line in table_parameters ] + [ [False]+line for line in table_parameters ]


table_parameters_mass = []
for mass in masses:
    table_parameters_mass = table_parameters_mass + [line+[mass] for line in table_parameters]

table_date = time.strftime("%d%m%y")
for line in table_parameters_mass:
    line[2] = table_date

def process_table_tab_for_pool(prams):
    process_table_tau.process_table_tau(*prams,taumax=7.)

def generate_h_weighted_table_for_pool(prams):
    h_weighted_table.generate_h_weighted_table(*prams)

with Pool(processes=64) as pool:
#     pool.map(process_table_tab_for_pool,table_parameters_taumode)
    pool.map(generate_h_weighted_table_for_pool,table_parameters_mass)

