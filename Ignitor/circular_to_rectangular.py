# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 10:33:36 2026

@author: felix
"""



import numpy as np

N_elem = 2
d_circ = 1.93

A = np.pi*d_circ**2/4

A_per_elem = A/N_elem

e = np.sqrt(A_per_elem)         # edge length per element

print(e)

e = 1.2

A = e**2*N_elem

d_circ = np.sqrt(4*A/np.pi)

print(d_circ)