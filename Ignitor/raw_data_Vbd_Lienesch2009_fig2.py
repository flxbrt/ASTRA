# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:19:58 2026

@author: felix
"""



import numpy as np

# concentration in % 
h2_conc = np.array([10.02194476,
                    14.98538021,
                    20.06486969,
                    25.02957117,
                    29.99776935,
                    35.02312038,
                    39.99348892,
                    45.13513514,
                    49.98932906,
                    55.01576526,
                    60.10550367,
                    65.02161318,
                    69.9930669,
                    75.07575164])

# breakdown voltage in kV converted to V (see below)
Vbd = np.array([14.49707475,
                12.95555841,
                11.6488631,
                10.65536299,
                11.17543151,
                11.43451839,
                12.89404331,
                13.49233527,
                14.66484712,
                15.3936622,
                18.52328875,
                21.49640383,
                23.42565694,
                23.50205022])*1e3

data = np.array([h2_conc, Vbd])

np.save('raw_data_Vbd_Lienesch2009_fig2.npy', data)