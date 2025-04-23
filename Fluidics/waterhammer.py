# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:06:11 2025

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi

# https://de.wikipedia.org/wiki/Drucksto%C3%9F

m_dot = 0.08
d = 0.01

l = 1 # longer than it will be

A = np.pi*d**2/4

p = 100e5
T = 293
fluid = 'O2'

rho = psi('D','P',p,'T',T,fluid)
sos = psi('A','P',p,'T',T,fluid)

print(sos, rho)

deltav = m_dot/rho/A

T_r = 2*l/sos

T_s = 10e-3 # faster than it will be

deltap = sos*rho*deltav*T_r/T_s/1e5

print(deltap) # overconservative due to longer pipe and faster valve than it actually will be 