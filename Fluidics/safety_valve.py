# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 18:39:03 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi

def get_q(Cv, p, T, fluid):
    # https://www.swagelok.com/downloads/webcatalogs/EN/MS-06-84.pdf
    # assuming chocked pressure reducer fully open
    # return volume flow in std l/min
    p = p/1e5
    
    rho = psi('D','P',1.013e5,'T',293,fluid)
    rho_air = psi('D','P',1.013e5,'T',293,'Air')
    Gg = rho/rho_air
    N2 = 6950
    
    return 0.471*N2*Cv*p/np.sqrt(Gg*T)

def calc_chocked_mass_flow(d_th, gamma, M, T, p):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        m_dot = np.pi*d_th**2/4*p/np.sqrt(T*R_ideal/M)*theta
        return m_dot

q_pressure_reducer = get_q(Cv=0.3, p=300e5, T=293, fluid='N2')
m_dot_pressure_reducer = psi('D','P',1.013e5,'T',293,'N2')*q_pressure_reducer/60/1000

# entweder gibt safety valve hersteller an, welchen volumenstrom sie bei welchem druck abführen können
# oder die geben öffnungsdurchmesser an und man muss sich den massenstrom ausrechnen

m_dot_safety_valve = calc_chocked_mass_flow(d_th=8e-3, gamma=1.4, M=28, T=270, p=40e5)

if m_dot_safety_valve > m_dot_pressure_reducer:
    print('Works :)')
else:
    print('no good')