# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 08:34:35 2026

@author: felix
"""



import numpy as np

# 1.0 atmosphere

# concentration in % 
h2_conc = np.array([12.32487876,
                    15.87819782,
                    20.44523277,
                    25.92159069,
                    32.88059927,
                    41.93798934,
                    54.24317721])

# quenching distance in mm converted to m (see below)
dq = np.array([1.951912568, 
                1.049180328, 
                0.681967213, 
                0.469945355, 
                0.469945355, 
                0.518032787, 
                0.843715847])*1e-3

data = np.array([h2_conc, dq])

np.save('raw_data_dq_Kim2001_fig3_1.0atm.npy', data)


# 1.5 atmosphere

# concentration in % 
h2_conc = np.array([12.22178812,
                    15.87543692,
                    20.34729928,
                    25.87189465,
                    32.92576328,
                    41.98320545,
                    54.19051191])

# quenching distance in mm converted to m (see below)
dq = np.array([1.626229508,
                0.933333333,
                0.572677596,
                0.384699454,
                0.365027322,
                0.415300546,
                0.633879781,
                ])*1e-3

data = np.array([h2_conc, dq])

np.save('raw_data_dq_Kim2001_fig3_1.5atm.npy', data)



# 2.0 atmosphere

# concentration in % 
h2_conc = np.array([12.26273266,
                    15.91992374,
                    20.39324467,
                    25.91685029,
                    32.92357541,
                    41.9814864,
                    54.28188178])

# quenching distance in mm converted to m (see below)
dq = np.array([1.344262295,
                0.8,
                0.500546448,
                0.271038251,
                0.273224044,
                0.343169399,
                0.467759563])*1e-3

data = np.array([h2_conc, dq])

np.save('raw_data_dq_Kim2001_fig3_2.0atm.npy', data)



# 2.5 atmosphere

# concentration in % 
h2_conc = np.array([12.25721087,
                    15.91705866,
                    20.39006704,
                    25.82001073,
                    32.92170008,
                    42.02638995,
                    54.28068366])

# quenching distance in mm converted to m (see below)
dq = np.array([1.112568306,
                0.679781421,
                0.367213115,
                0.207650273,
                0.194535519,
                0.227322404,
                0.417486339])*1e-3

data = np.array([h2_conc, dq])

np.save('raw_data_dq_Kim2001_fig3_2.5atm.npy', data)