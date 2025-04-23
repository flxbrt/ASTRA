# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 15:59:47 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi

# calculation for pressurizer downstream the reducer

m_dot = 1.1464120324505944 # max mass flow at full thrust
rof = 1.1
m_dot_ox = rof/(1+rof)*m_dot
m_dot_f = 1/(1+rof)*m_dot

p_tank = 40e5
T_ox = 120
T_f = 300

rho_ox = psi('D','P',p_tank,'T',T_ox,'O2')
rho_f = psi('D','P',p_tank,'T',T_f,'C2H6O')

vol_dot_ox = m_dot_ox/rho_ox
vol_dot_f = m_dot_f/rho_f

# single lines for ox and fuel

T_press = 270 ### highest temp results in lowest density --> after reducer, so lower than in tank due to joule thompson
p_press = 40e5

rho_press = psi('D','P',p_press,'T',T_press,'N2')

d_pipe_ox = 4e-3
d_pipe_f = 4e-3

A_ox = d_pipe_ox**2/4*np.pi
A_f = d_pipe_f**2/4*np.pi

vol_dot_press_ox = vol_dot_ox
vol_dot_press_f = vol_dot_f

vel_press_ox = vol_dot_press_ox/A_ox
vel_press_f = vol_dot_press_f/A_f

print(f'{vel_press_ox=} [m/s] \n{vel_press_f=} [m/s]')

speed_of_sound = psi('speed_of_sound','P',p_press,'T',T_press,'N2')

print(f'{speed_of_sound=} [m/s]')

p_dyn_ox_approx = 0.5*rho_press*vel_press_ox**2/1e5
p_dyn_f_approx = 0.5*rho_press*vel_press_f**2/1e5

print(f'{p_dyn_ox_approx=} [bar] \n{p_dyn_f_approx=} [bar]')

m_dot_press = (vol_dot_press_ox+vol_dot_press_f)*rho_press

# calculation for pressurizer upstream the reducer

# more critical is the condition, where the density is lowest --> highest velocity

T_press = 220
p_press = 100e5

rho_press_upstream =  psi('D','P',p_press,'T',T_press,'N2')

vol_dot_press_upstream = m_dot_press / rho_press_upstream

d_pipe = 4e-3

A_pipe = d_pipe**2/4*np.pi

vel_upstream = vol_dot_press_upstream/A_pipe

print(f'{vel_upstream=} [m/s]')

speed_of_sound = psi('speed_of_sound','P',p_press,'T',T_press,'N2')

print(f'{speed_of_sound=} [m/s]')