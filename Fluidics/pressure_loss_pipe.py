# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 12:21:27 2024

@author: felix
"""



import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as psi


def reynolds(vel, diameter, rho, eta):
    
    return vel*diameter*rho/eta

def friction_factor(Re, k=0):
    if Re <= 2320:
        lam = 64/Re
        print(1)
    elif Re > 2320 and Re <= 1e5:
        lam = 0.316/(Re**(0.25))
        print(2)
    elif Re > 1e5:
    # else:
        lam_up = 10
        lam_low = 0
        err = 100
        while abs(err) > 1e-4:
            lam_mid = np.mean([lam_low, lam_up])
            err = 1/np.sqrt(lam_mid)+2*np.log10(2.51/Re/np.sqrt(lam_mid) + k/3.71/d)
            if err>0:
                lam_low = lam_mid
            elif err<0:
                lam_up = lam_mid
            elif err == 0:
                break
        lam = lam_mid
        print(3)
    return lam


# https://www.schweizer-fn.de/stroemung/rauhigkeit/rauhigkeit.php 

#%%

dp_ethanol = list()
dp_lox = list()
l = 1

for d in np.linspace(8e-3, 10e-3, 50):
    p = 40e5
    T = 293
    
    rho = psi('D','P',p,'T',T,'C2H6O')
    eta = psi('viscosity','P',p,'T',T,'C2H6O')

    
    m_dot = 0.73
    
    A = np.pi*d**2/4
    
    vel = m_dot/A/rho
    
    Re = reynolds(vel, d, rho, eta)
    
    lam = friction_factor(Re, 50e-6)
    
    deltap_ethanol = lam*l*rho*vel**2/(2*d)/1e5
    
    p = 40e5
    T = 110
    rho = psi('D','P',p,'T',T,'O2')
    eta = psi('viscosity','P',p,'T',T,'O2')
    
    m_dot = 0.85
    
    A = np.pi*d**2/4
    
    vel = m_dot/A/rho
    
    Re = reynolds(vel, d, rho, eta)
    
    lam = friction_factor(Re, 50e-6)
    
    deltap_lox = lam*l*rho*vel**2/(2*d)/1e5
    
    dp_ethanol.append(deltap_ethanol)
    dp_lox.append(deltap_lox)

plt.plot(np.linspace(8, 10, 50), dp_ethanol, label='Ethanol')
plt.plot(np.linspace(8, 10, 50), dp_lox, label='Oxygen')
plt.legend()
plt.xlabel('Pipe Diameter [mm]')
plt.ylabel('Pressure Loos PER METER [bar/m]')
plt.grid()



#%%

# cooling schannel estimation

h = 1e-3
w = 1e-3/2
A = h*w
U = 2*(h+w)

d_hyd = 4*A/U



m_dot = 0.73/40 # 0.687/52 



# d_hyd = 0.001
l = 0.2
# m_dot = 10
A = d_hyd**2/4*np.pi

rho = psi('D','P',40e5,'T',293,'C2H6O')
eta = psi('viscosity','P',40e5,'T',293,'C2H6O')

vel = m_dot/A/rho

# l = 0.3



Re = reynolds(vel, d_hyd, rho, eta)

lam = friction_factor(Re, k=50e-6)

# lam = 0.0715

deltap = lam*l*rho*vel**2/(2*d_hyd)/1e5

# deltap2 = lam*l/d/2/A**2/rho*m_dot**2/1e5

#%% hihg pressure system


dp_ethanol = list()
dp_lox = list()
l = 1

for d in np.linspace(6e-3, 10e-3, 50):
    p = 35e5
    T = 293
    
    fluid = 'N2'
    
    rho = psi('D','P',p,'T',T,fluid)
    eta = psi('viscosity','P',p,'T',T,fluid)
    
    Q_std = 4900 # l/min
    Q_std = Q_std/1000/60 # m^3/s
    rho_std = psi('D','P',1e5,'T',293,fluid)
    
    m_dot = rho_std*Q_std
    
    A = np.pi*d**2/4
    
    vel = m_dot/A/rho
    
    Re = reynolds(vel, d, rho, eta)
    
    lam = friction_factor(Re, 50e-6)
    
    deltap_ethanol = lam*l*rho*vel**2/(2*d)/1e5
    
    Q_std = 5300 # l/min
    Q_std = Q_std/1000/60 # m^3/s
    rho_std = psi('D','P',1e5,'T',293,fluid)
    
    m_dot = rho_std*Q_std
    
    vel = m_dot/A/rho
    
    Re = reynolds(vel, d, rho, eta)
    
    lam = friction_factor(Re, 50e-6)
    
    deltap_lox = lam*l*rho*vel**2/(2*d)/1e5
    
    dp_ethanol.append(deltap_ethanol)
    dp_lox.append(deltap_lox)

plt.plot(np.linspace(8, 10, 50), dp_ethanol, label='Ethanol')
plt.plot(np.linspace(8, 10, 50), dp_lox, label='Oxygen')
plt.legend()
plt.xlabel('Pipe Diameter [mm]')
plt.ylabel('Pressure Loos PER METER [bar/m]')
plt.grid()

#%% nitrogen

fluid = 'N2'

d = 6e-3

Q_std = 6200 # std l / min

rhoN = psi('D','P',1e5,'T',293,fluid)

m_dot = Q_std/60/1000*rhoN

l = 4

A = d**2/4*np.pi

p = 150e5
T = 293

rho = psi('D','P',p,'T',T,fluid)
eta = psi('viscosity','P',p,'T',T,fluid) # dynamic viscosity
nu = eta/rho

vel = m_dot/A/rho

Re = reynolds(vel, d, rho, eta)

lam = friction_factor(Re, k=50e-6)

deltap = lam*l*rho*vel**2/(2*d)/1e5

#%% injector


# general 
ROF = 1.17
mdot = 1.26

mdot_fuel = 1/(1+ROF) * mdot
mdot_ox = ROF/(1+ROF) * mdot

# LOx
fluid = 'O2'
d = 4e-3
N = 8

mdot_ox = mdot_ox/N

l = 0.01

A = d**2/4*np.pi

p = 30e5
T = 90

rho = psi('D','P',p,'T',T,fluid)
eta = psi('viscosity','P',p,'T',T,fluid) # dynamic viscosity
nu = eta/rho

vel = m_dot/A/rho

Re = reynolds(vel, d, rho, eta)

lam = friction_factor(Re, k=50e-6)

deltap1 = lam*l*rho*vel**2/(2*d)/1e5
deltap2 = 0.5*rho*vel**2/1e5
# Ethanol
