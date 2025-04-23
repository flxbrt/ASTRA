# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 19:18:59 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi



# WTP Skript SS22 p. 50f.

#%% Tank

V = 0.1
d = 0.3
A = d**2/4*np.pi
h = V/A

A_wetted = d*np.pi*h

p_O2 = 50e5

Tw = T_lox = psi('T','P',p_O2,'Q',0,'O2')

# worst case: LOx at 90K
Tw = 90

L = h

#%% Rohrleitung LOx Befüllung

l = 3
d = 0.01
A = d**2/4*np.pi
A_wetted = d*np.pi*l

Kf = 0.8
Nu0 = 0.36

#%%

# worst case: LOx at 90K
Tw = 90

L = d


T_amb = 303
beta = 1/T_amb
g = 9.81
# Tw = 90




rho = psi('D','P',1e5,'T',T_amb,'Air')
eta = psi('viscosity','P',1e5,'T',T_amb,'Air')
lam = psi('conductivity','P',1e5,'T',T_amb,'Air')



nu = eta/rho


Pr = 0.7                                # https://www.chemie-schule.de/KnowHow/Prandtl-Zahl
f_pr = 0.913
Gr = g*beta*L**3*(Tw-T_amb)/nu**2


Ra = abs(Gr*Pr)

if Ra > 1e9:
    Nu = 0.15*f_pr**(4/3)*Ra**(1/3)
elif Ra > 1e4 and Ra < 1e9:
    Nu = Nu0 + 0.668*Kf*f_pr*Ra**0.25
else:
    print('No correlation is valid')
    
alpha = Nu * lam / L


Q = alpha*A_wetted*(T_amb-Tw)

# sehr konservative abschätzung da angenommen wird, dass wand und fluid an der innenseite
# keinen wärmewiderstand hat 
# das führt dazu, dass angenommen wird, dass an der außenseite die LOx temperatur herrscht

h1 = psi('H','P',p_O2,'Q',0,'O2') # J/kg
h2 = psi('H','P',p_O2,'Q',1,'O2')

m_dot = Q/(h2-h1)

print(m_dot)

rhoN = psi('D','P',1e5,'T',293,'O2')

V_dot = m_dot/rhoN

print(V_dot*1000/60)

#%%

p_O2 = 1e5

p_open = 96e5

deltat = 0.1

rhoL = psi('D','P',p_O2,'Q',0,'O2')

m = V*rhoL

dh_dt = Q/m

h1 = psi('H','P',p_open,'D',rhoL,'O2')

h2 = h1+dh_dt*deltat

drho_dt = (psi('D','P',p_open,'H',h1,'O2')-psi('D','P',p_open,'H',h2,'O2'))/deltat

m_dot = drho_dt*V