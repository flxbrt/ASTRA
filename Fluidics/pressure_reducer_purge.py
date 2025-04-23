# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:10:44 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi



def calc_chocked_mass_flow(A, gamma, M, T, p):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        m_dot = A*p/np.sqrt(T*R_ideal/M)*theta
        return m_dot
    
    
A_ox = 1.45/1.15*21.991e-6

A_fu = 1.45/1.15*20.232e-6

gamma = 1.4
M = 28
T = 293
p = 30e5

# A = A_ox


d = 9.5e-3
A = d**2/4*np.pi

# A = 8e-3**2/4*np.pi

m_dot = calc_chocked_mass_flow(A, gamma, M, T, p)

rhoN = psi('D','P',1e5,'T',293,'N2')

q_std = m_dot/rhoN
q_std_l_min = q_std*1000*60
print(m_dot, q_std_l_min)

#%%

fluid = 'O2'
fluid = 'C2H6O'

dil_fac = 1

if fluid == 'O2':
    m_dot_l = 0.85
    p_l = 30e5 + 9e5
    T_l = 110
if fluid == 'C2H6O':
    m_dot_l = 0.73
    p_l = 30e5 + 9e5 + 9e5
    T_l = 288

rho_l = psi('D','P',p_l,'T',T_l,fluid)

Q_l = m_dot_l/rho_l

Q_p = Q_l*dil_fac

rho_p = psi('D','P',p_l,'T',293,'N2')

m_dot_p = Q_p * rho_p

rhoN = psi('D','P',1e5,'T',293,'N2')

Q_std = m_dot_p / rhoN * 1000*60