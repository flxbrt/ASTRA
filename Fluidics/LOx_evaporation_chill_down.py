# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 12:03:44 2025

@author: felix
"""



from CoolProp.CoolProp import PropsSI as psi

cp_steel = 477
m_steel = 50

T_lox = 90
T_amb = 293

deltaT = T_amb - T_lox

h_li = psi('H','P',1e5,'Q',0,'O2')
h_ga = psi('H','P',1e5,'Q',1,'O2')

deltaH = h_ga-h_li

E_steel = m_steel*cp_steel*deltaT

# E_lox = m_lox*deltaH

# V_lox = 10*1e-3
# m_lox = rho_lox*V_lox

m_lox = E_steel/deltaH

rho_lox = psi('D','P',1e5,'Q',0,'O2')

V_lox = m_lox/rho_lox*1000