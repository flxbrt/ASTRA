# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 16:51:07 2024

@author: felix -> felix.ebert@tum.de

This script implements a one-dimensional radial transient heat transfer 
problem in a hollow cylinder in which a fuel and oxidiser react and combust
It is intended to be used for a parametric study on the time scale until 
the wall temperature threshold is reached for various setups (geometry, 
wall material, propellant combination etc.)

Main assumptions
- geometry: hollow cylinder
- 1D radial heat transfer only
- inner side heat transfer via convection only
- no outside heat transfer, i.e. all the heat that goes into the wall stays
in the wall material
- homogenous flow in the cross sectinal plane of the hollow cylinder

temporal discretisation via RK4
spatial discretisation via central difference except for the boundary conditions
"""

import numpy as np
import matplotlib.pyplot as plt



# https://pages.mtu.edu/~fmorriso/cm310/energy2013
# conservation of energy in differential form neglecting (1) convective heat
# transport and (2) circumenferentiel and axial heat transport due to conduction
# \rho c_p \frac{\partial T}{\partial t} = -\frac{1}{r} 
        # \frac{\partial}{\partial r} \left( r q \right)





def rk4(f, y, h, t, *args):
    # runge kutte 4th order explicit
    tk_05 = t + 0.5*h
    yk_025 = y + 0.5 * h * f(t, y, *args)
    yk_05 = y + 0.5 * h * f(tk_05, yk_025, *args)
    yk_075 = y + h * f(tk_05, yk_05, *args)
    
    return y + h/6 * (f(t, y, *args) + 2 * f(tk_05, yk_025, *args) + 2 * f(tk_05, yk_05, *args) + f(t+h, yk_075, *args))





def time_derivative(t, y, *args):
        
    radial = args[0]
    dr = args[1]
    rho = args[2]
    cp = args[3]
    lam = args[4]
    q_hg = args[5]
    
    # previous solution
    Temperature = y
    
    # next = outflow
    # T_next = np.roll(Temperature, shift=-1)
    # prev = inflow
    # T_prev = np.roll(Temperature, shift=1)
    
    q = np.zeros(len(Temperature)+1)
    q[1:-1] = -lam * (Temperature[1:] - Temperature[:-1])/dr
    q[0] = q_hg
    q[-1] = 0
    
    q_in = q[:-1]
    q_out = q[1:]
    # q = -lam * (T_next - T_prev)/dr
    # q_next = np.roll(q, shift=-1)
    # q_prev = np.roll(q, shift=-1)
    
    # next = outflow
    r_next = np.roll(radial, shift=-1)
    # prev = inflow
    r_prev = np.roll(radial, shift=1)
    
    dT_dt = 1/radial/rho/cp*((r_prev*q_in - r_next*q_out)/2/dr)
    
    # boundary conditions
    # dT_dt[0] = 1/radial[0]/rho/cp*((r_prev[0]*q_hg - r_next[0]*q_next[0])/2/dr)
    # dT_dt[-1] = 1/radial[-1]/rho/cp*((r_prev[-1]*q_prev[-1] - 0)/2/dr)
    
    # print(q)
    
    return np.array(dT_dt)
    



# x integrator
# x dynamik
# plotting
# verbrennung vorgeben
# x geometrie und materialdaten vorgeben





ich glaub hier liegt im wesentlichen ein problem bei der räumlichen diskreitisuerng vor
entweder ich arbeite mich nochmal durh die diskretisierungsschemate durch oder
ich stelle eine energiebeilanz inkl. der flächen auf und rechne damit --> also eine integrale energiebilanz





# https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=NINC34
rho = 8190
cp = 435
lam = 11.4

r_inner = 0.01
r_outer = 0.02
dr = 1e-5
radial = np.arange(r_inner, r_outer, dr)

T_init = 300

hg = 1e5
T_aw = 2000

t_start = 0
t_end = 1
h = 1e-5

time = np.arange(t_start, t_end, h)

temperature_array = np.ones((len(time), len(radial)))*T_init

y = temperature_array[0,:]


#%%
if __name__ == '__main__':
    
    for t_index, t in enumerate(time):
        
        q_hg = hg*(T_aw-y[0])
        
        args = [radial, dr, rho, cp, lam, q_hg]
        
        y = rk4(time_derivative, y, h, t, *args)
        
        temperature_array[t_index, :] = y


#%%
plt.plot(radial, temperature_array[0,:])
plt.plot(radial, temperature_array[50,:])
plt.plot(radial, temperature_array[99,:])
plt.xlabel('Radius [m]')
plt.ylabel('Temperature [K]')


#%%

print(lam * h / (rho * cp * dr**2))
