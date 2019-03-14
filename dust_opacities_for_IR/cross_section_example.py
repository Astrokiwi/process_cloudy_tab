import numpy as np
import cloud_cross_section

cc_func = cloud_cross_section.get_cloud_cross_section_function(file="data/constant_density_table_120319.dat")

bb_temp = np.log10(150.)
NH_cent = 10.**25 # cm**-2

surf = cloud_cross_section.full_NH_to_table_density(NH_cent)

cross_section = cc_func(surf,bb_temp)

print("result=",cross_section)

# say, cloud is 0.01 pc in radius

cloud_rad = 0.01 * cloud_cross_section.pc
cloud_area = np.pi * cloud_rad**2

print("Cloud geometric area:",cloud_area," cm**2")
print("Cloud cross section:",cloud_area*(10.**cross_section)," cm**2")