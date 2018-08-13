import process_table_tau
import h_weighted_table
import time

for nodustmode,highdensemode,filename in [ (False,True,"highden_260118.txt"),
                                           (False,False,"tables_070818.txt"),
                                           (True,False,"nodust_301117.txt")
                                            ]:
    for taumode in [True,False]:
        process_table_tau.process_table_tau(taumode,nodustmode,highdensemode,filename,taumax=7.)
    table_date = time.strftime("%d%m%y")
#     print(table_date)
    h_weighted_table.generate_h_weighted_table(nodustmode,highdensemode,table_date,mass=0.0001)
