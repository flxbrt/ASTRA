# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 16:48:48 2024

This script implements a one-dimensional transient radial heat transfer 
problem in a hollow cylinder in which a fuel and oxidiser react and combust

Main assumptions
- geometry: hollow cylinder
- 1D radial heat transfer
--> symmetry in circumferential direction
--> no axial heat flow
- inner side heat transfer via convection only --> no radiative heat transfer
- no outside heat transfer --> adiabativ boundary condition at the outer wall
of the hollow cylinder
- homogenous flow in the cross sectinal plane of the hollow cylinder

temporal discretisation via RK4
spatial discretisation via finite / central difference

@author: felix
"""

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

def rk4(f, y, h, t, *args):
    # runge kutte 4th order explicit
    tk_05 = t + 0.5*h
    yk_025 = y + 0.5 * h * f(t, y, *args)
    yk_05 = y + 0.5 * h * f(tk_05, yk_025, *args)
    yk_075 = y + h * f(tk_05, yk_05, *args)
    
    return y + h/6 * (f(t, y, *args) + 2 * f(tk_05, yk_025, *args) + 2 * f(tk_05, yk_05, *args) + f(t+h, yk_075, *args))

def set_combustion(p, fuel, ox, ROF):
    comb = ct.Solution('gri30_WARR.yaml')
    # comb = ct.Solution('gri30_WARR.yaml')
    comb.Y = {fuel: 1, ox: ROF}
    comb.TP = 273, p
    comb.equilibrate('HP')
    return comb

def get_hg(comb, d, A, p, T_w, m_dot):
    ''' Bartz equation implemented '''
    T0 = comb.T
    gamma = comb.cp / comb.cv
    a = np.sqrt(gamma * ct.gas_constant * comb.T / comb.mean_molecular_weight)
    mu = comb.viscosity 
    lam = comb.thermal_conductivity
    Pr = (mu * cp) / lam
    c_star = p*d**2/4*np.pi/m_dot
    vel = m_dot / A / comb.density
    # print(vel)
    Ma = vel / a
    # print(Ma)
    C = 0.026
    omega = 0.6
    sigma = 1 / (((((1 / 2) * (T_w / T0) * (1 + (((gamma - 1) / 2) * Ma**2))) + (1 / 2))**(0.8 - (omega / 5))) * ((1 + (((gamma - 1) / 2) * Ma**2))**(omega / 5)))
    hg = (C / (d**0.2)) * ((mu**0.2 * cp)/(Pr**0.6)) * ((p / c_star)**0.8) * sigma # * fac
    
    return hg

def time_derivative(t, y, *args):
    ''' 
    differential equations for conservation of energy 
    dU/dt = dQ/dt
    dQ/dt --> obtained per node through convective heat flux or through heat conduction
    '''
    dr = args[0]
    cp = args[1]
    lam = args[2]
    m = args[3]
    Q_i = args[4]
    Q_e = args[5]
    A = args[6]
    
    # solution rfrom previous time step
    Q = y[0,:]
    T = y[1,:]
    
    dQ_dt = np.zeros(len(Q))
    # BC on hot gas side
    dQ_dt[0] = Q_i - (-lam*(T[1] - T[0])/dr*A[1])
    # heat conduction
    dQ_dt[1:-1] = -lam/dr*(T[1:-1] - T[:-2])*A[1:-2] - (-lam/dr*(T[2:] - T[1:-1])*A[2:-1])
    # BC on adiabatic side
    dQ_dt[-1] = -lam*(T[-1] - T[-2])/dr*A[-2] - Q_e
    
    dT_dt = dQ_dt/m/cp
    
    return np.array([dQ_dt, dT_dt])   

#%% problem definition

# Material properties: 1.4571
rho = 8000
cp = 500
lam = 15
T_melt = 1398 + 273
T_usable = 700 + 273
# Material properties: CCZr
# rho = 8900
# cp = 394
# lam = 335
# T_melt = 1070 + 273
# T_usable = 500 + 273

# Geometry
l_unit = 1          # unit length of cylinder --> can be set to arbitrary value --> only required as senity checks to match physical units
r_inner = 1      # inner radius of hollow cylinder
r_outer = 1.02      # outer radius of hollow cylinder
dr = 1e-4           # spatial discretisation
radial_c = np.arange(r_inner, r_outer + 1e-6, dr) # radius at respective node centre
radial_e = np.arange(r_inner-dr/2, r_outer+dr/2 + 1e-6, dr) # radius at node edge / interface between two nodes

A = 2*np.pi*radial_e*l_unit     # cross section perpendicular to heat flux at each node edge
V = (radial_e[1:]**2-radial_e[:-1]**2)*np.pi*l_unit # volume of each node / element
m = rho*V                       # mass of each element

# Propellant
# fuel = 'C2H5OH'
# ox = 'O2'
# ROF = 1.4                     # [-]
# p = 20e5                    # [Pa]
# m_dot = 1.2              # [kg]

# comb = set_combustion(p, fuel, ox, ROF)

# Initial Conditions and Integration Parameters
T_init = 300                # Initial Temperature
t_start = 0                 # Simulation start time
t_end = 6.5                   # Simulation end time
h = 1e-3                    # time step size

time = np.arange(t_start, t_end, h)

# Solution arrays
temperature_array = np.ones((len(time)+1, len(radial_c)))*T_init
internal_heat_array = np.ones((len(time)+1, len(radial_c)))*T_init*cp*m

#%% Main Loop

y = np.array([internal_heat_array[0,:], temperature_array[0,:]])

# T_aw = comb.T       # actually adiabatic wall temperature lower than stagnation
#                     # temperature --> this assumption only holds for low flow velocity
    
# throat
# hg = 11254                
# T_aw = 2767
# faceplate
hg = 933
T_aw = 2826

if __name__ == '__main__':
    for t_index, t in enumerate(time):
        T_w = y[1,0]
        # hg = get_hg(comb, r_inner*2, r_inner**2*np.pi, p, T_w, m_dot)         # convective heat transfer coefficient
        # print(hg)
        # the get_hg functions shouldn't be used
        # it utilizes too many simplifications that are not valid
        
        Q_i = hg*(T_aw-T_w)*A[0]                    # inner BC heat flux into the wall
        Q_e = 0                                     # outer BC heat flux exiting the wall
        
        # print(hg*(T_aw-T_w)/1e6)
        # if t_index == 1:
        #     break
        
        args = [dr, cp, lam, m, Q_i, Q_e, A]
        y = rk4(time_derivative, y, h, t, *args)    # time integration of all radial nodes at once
        
        internal_heat_array[t_index+1, :] = y[0,:]
        temperature_array[t_index+1, :] = y[1,:]

#%% Numerical Stability
# condition for numerical stability
# should be a bit lower than 0.5 to obtain numerically stable result

condition = lam * h / (rho * cp * dr**2)
print(f'Simulation is numerically stable: {condition<0.5}')
#%% Plotting
# plot temperature in wall at last time step

plt.plot(radial_c, temperature_array[-1,:])
plt.hlines(y=T_melt, xmin=r_inner, xmax=r_outer, color='orange', label="Melting Point")
plt.hlines(y=T_usable, xmin=r_inner, xmax=r_outer, label="Usable")
plt.legend()
plt.grid()
plt.xlabel('Radius [m]')
plt.ylabel('Temperature [K]')