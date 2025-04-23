# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 14:05:18 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi

pressure = np.linspace(10e5, 40e5, 100)

T_sat = list()
rho_sat = list()

for p in pressure:
    T_sat.append(psi('T','P',p,'Q',1,'O2'))
    rho_sat.append(psi('D','P',p,'Q',0,'O2'))
    
    
# --> 700 kg/m^3 is lowest density we might have --> results in highest LOx velocity in pipe

m_dot = 1.1464120324505944
ROF = 1.1
m_dot_lox = ROF/(1+ROF)*m_dot

rho = 850

diameter = np.linspace(4e-3, 8e-3, 100)

vel = list()

for d in diameter:
    A = d**2/4*np.pi
    vel.append(m_dot_lox/rho/A)    
    
import matplotlib.pyplot as plt

plt.plot(diameter*1000, vel)