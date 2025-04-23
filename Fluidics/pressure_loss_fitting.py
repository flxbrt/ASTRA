# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:02:14 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi


# d2 = 9.5mm --> tee or cross inner diameter
# d1 = 10mm --> pipe
# d2/d1 = 0.95 --> zeta = 0.05 https://www.schweizer-fn.de/zeta/verengung/verengung.php

zeta = 0.4
d = 0.01
l = 0.03

# 90° bend with R/D = 1
# https://www.schweizer-fn.de/zeta/rohrbogen/rohrbogen.php
# https://www.schweizer-fn.de/berechnung/zeta/rohrbogen/rohrbogen_rech.php

# zeta = 0.225
# d = 0.01
# l = 0.05

p = 40e5
T = 100
rho = psi('D','P',p,'T',T,'O2')

m_dot = 0.8

A = np.pi*d**2/4

vel = m_dot/A/rho

# Armaturen, Formstücke https://www.schweizer-fn.de/stroemung/druckverlust/druckverlust.php#druckverlustzeta

deltap = zeta*rho*vel**2/2/1e5 * l/d



#%% Aufweitung 300bar Schlauch --> von 3.5 auf 6mm
zeta = 0.6 # konservativ
d = 3.5e-3

rho = psi('D','P',1.013e5,'T',293,'N2')
V_dot = 6200/1000/60

m_dot = rho*V_dot

# m_dot = 1

p = 150e5
T = 293
rho = psi('D','P',p,'T',T,'N2')

A = np.pi*d**2/4

vel = m_dot/A/rho

deltap = zeta*rho*vel**2/2/1e5 #* l/d

def calc_chocked_mass_flow(d_th, gamma, M, T, p):
    # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
    theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
    R_ideal = 8314
    m_dot = np.pi*d_th**2/4*p/np.sqrt(T*R_ideal/M)*theta
    return m_dot

m_dot_chocked = calc_chocked_mass_flow(d, 1.4, 28, T, p)

eta = psi('viscosity','P',p,'T',T,'N2')
