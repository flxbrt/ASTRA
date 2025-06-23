# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:31:59 2025

@author: felix
"""


import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as psi


# rausfinden, wie zeta mit massenstrom skaliert


def get_fluid_properties(p, T, fluid):
    rho = psi('D', 'P', p, 'T', T, fluid)
    sos = psi('speed_of_sound', 'P', p, 'T', T, fluid)
    eta = psi('V', 'P', p, 'T', T, fluid)
    return rho, sos, eta


def friction_factor(massflow, rho, cross_section, length, diameter, eta):
    vel = massflow/rho/cross_section
    Re = vel*diameter*rho/eta
    # CM Diss p.27
    # if self.Re == 0:
    # self.lam = 10000
    if abs(Re) <= 2320:
        if abs(Re) < 1:
            lam = 64
        else:
            lam = 64/abs(Re)
        # if self.lam > 1500:
            # self.lam = 1500
        # print(1)
        corr = 1
    elif abs(Re) > 2320 and abs(Re) <= 1e5:
        lam = 0.316/(abs(Re))**(0.25)
        # print(2)
        corr = 2
    elif abs(Re) > 1e5:
        lam_up = 10
        lam_low = 0
        err = 100
        while abs(err) > 1e-4:
            lam_mid = np.mean([lam_low, lam_up])
            err = 1/np.sqrt(lam_mid)+2*np.log10(2.51/abs(Re)/np.sqrt(lam_mid))
            if err > 0:
                lam_low = lam_mid
            elif err < 0:
                lam_up = lam_mid
            elif err == 0:
                break
        lam = lam_mid
        # print(3)
        corr = 3
    zeta = lam*length/diameter/(2*cross_section**2)
    return zeta, corr


d = 0.001
l = 1
cross_section = np.pi*d**2/4

p_f = 40e5
T_f = 293
fluid = 'C2H6O'

rho, sos, eta = get_fluid_properties(p_f, T_f, fluid)
massflow = np.linspace(0.3, 0.7, 100)

zeta = np.zeros(len(massflow))
correlation = np.zeros(len(massflow))

for index, m_dot in enumerate(massflow):
    zeta[index], correlation[index] = friction_factor(
        m_dot, rho, cross_section, l, d, eta)


try:
    one = np.where(correlation==1)[0][-1]
except IndexError:
    one=0
try:
    two = np.where(correlation==2)[0][-1]
except IndexError:
    two=0

# one = np.where(correlation==1)[0][-1]
# two = np.where(correlation==2)[0][-1]
thr = np.where(correlation==3)[0][-1]

plt.plot(massflow[:one], zeta[:one], color='b')
plt.plot(massflow[one:two], zeta[one:two], color='orange')
plt.plot(massflow[two:thr], zeta[two:thr], color='r')


#%% corr 2 und corr 3 übereinander plotten

def friction_factor_2(massflow, rho, cross_section, length, diameter, eta):
    vel = massflow/rho/cross_section
    Re = vel*diameter*rho/eta
    # CM Diss p.27
    # if self.Re == 0:
    # self.lam = 10000
    
    lam_2 = 0.316/(abs(Re))**(0.25)
    
    lam_up = 10
    lam_low = 0
    err = 100
    while abs(err) > 1e-4:
        lam_mid = np.mean([lam_low, lam_up])
        err = 1/np.sqrt(lam_mid)+2*np.log10(2.51/abs(Re)/np.sqrt(lam_mid))
        if err > 0:
            lam_low = lam_mid
        elif err < 0:
            lam_up = lam_mid
        elif err == 0:
            break
    lam_3 = lam_mid
    
    zeta_2 = lam_2*length/diameter/(2*cross_section**2)
    zeta_3 = lam_3*length/diameter/(2*cross_section**2)
    return zeta_2, zeta_3

zeta = np.zeros((2,len(massflow)))


for index, m_dot in enumerate(massflow):
    zeta[:,index] = friction_factor_2(
        m_dot, rho, cross_section, l, d, eta)
    
plt.plot(massflow, zeta[0,:], label='2')
plt.plot(massflow, zeta[1,:], label='3')
plt.legend()
