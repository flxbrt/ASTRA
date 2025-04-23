# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 17:20:46 2024

@author: felix
"""



import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as psi




# INTERPRETATION ERGEBNISSE
# WENN MAN AM ENDE DES SCHUBEVENTS NOCH DEN VOLLEN SCHUB ZUR VERFÜGUNG HABEN MÖCHTE, MUSS MAN MIT MAXIMALEM MASSENSTROM RECHNEN




def get_cv(q, p, T, fluid):
    # https://www.swagelok.com/downloads/webcatalogs/EN/MS-06-84.pdf
    
    p = p/1e5
    
    rho = psi('D','P',1.013e5,'T',293,fluid)
    rho_air = psi('D','P',1.013e5,'T',293,'Air')
    Gg = rho/rho_air
    N2 = 6950
    return q/0.471/N2/p*np.sqrt(Gg*T)


cv_ox = list()
cv_fuel = list()
q_std_fuel_list = list()
q_std_ox_list = list()

# hopper
# vol_press = 0.04107810793497826
# vol_fuel = 0.05613619288644423
# vol_ox = 0.043832435999593126
# m_dot_prop = 1.3268113449977343
# ROF = 1.025
# m_dot_fuel = 1/(1+ROF)*m_dot_prop
# m_dot_ox = ROF/(1+ROF)*m_dot_prop
# p_reducer_fuel = 40e5
# p_reducer_ox = 30e5
# flag = 'adiabatic'

# test bench
vol_press = 12*0.05
vol_fuel = 0.4
vol_ox = 0.4
m_dot_fuel = 0.73
m_dot_ox = 0.85
p_reducer_fuel = 87e5
p_reducer_ox = 79e5
flag = 'adiabatic' # 'isothermal'



ullage_fuel = 0.1
ullage_ox = 0.1

vol_ullage_fuel = vol_fuel*ullage_fuel
vol_ullage_ox = vol_ox*ullage_ox

vol_fuel_tank = vol_fuel + vol_ullage_fuel
vol_ox_tank = vol_ox + vol_ullage_ox

press = 'N2'

T_init = 273
p_init = 300e5
# p_end = 150e5




rho_press_init = psi('D','P',p_init,'T',T_init,press)
mass_press = rho_press_init*vol_press


# m_dot_prop = 1.2
# ROF = 1.1
# m_dot_fuel = 1/(1+ROF)*m_dot_prop
# m_dot_ox = ROF/(1+ROF)*m_dot_prop




T_fuel = 293
T_ox = 100
fuel = 'C2H6O'
ox = 'O2'
rho_fuel = psi('D','P',p_reducer_fuel,'T',T_fuel,fuel)
rho_ox = psi('D','P',p_reducer_ox,'T',T_ox,ox)

kappa = 1.4

# T_end = T_init*(p_end/p_init)**((kappa-1)/kappa)

# calculation
vol_dot_fuel = m_dot_fuel / rho_fuel
vol_dot_ox = m_dot_ox / rho_ox

N = 1000
deltat = vol_fuel / vol_dot_fuel # same result for deltat obtained when vol_ox/vol_dot_ox is used
dt = deltat / N

p_current = p_init

mass_press_array_jt = np.zeros(N)
T_reducer_fuel_array_jt = np.zeros(N)
T_reducer_ox_array_jt = np.zeros(N)
T_presstank_array_jt = np.zeros(N)

# rho_reducer_array = np.zeros(N)
p_press_array_jt = np.zeros(N)

vol_dot_press_fuel = vol_dot_fuel
vol_dot_press_ox = vol_dot_ox

for ii in range(N):
    if flag == 'adiabatic':
        T_current = T_init*(p_current/p_init)**((kappa-1)/kappa)
    elif flag == 'isothermal':
        T_current = T_init
    h_current = psi('H','P',p_current,'T',T_current,press)
    # fuel
    T_reducer_fuel = psi('T','P',p_reducer_fuel,'H',h_current,press)
    rho_reducer_fuel = psi('D','P',p_reducer_fuel,'H',h_current,press)
    # ox
    T_reducer_ox = psi('T','P',p_reducer_ox,'H',h_current,press)
    rho_reducer_ox = psi('D','P',p_reducer_ox,'H',h_current,press)
    
    m_dot_reducer = vol_dot_press_fuel*rho_reducer_fuel + vol_dot_press_ox*rho_reducer_ox
    
    
    '''
    
    
    
    
    
    
    
    
    
    
    '''
    m_dot_reducer = vol_dot_press_ox*rho_reducer_ox
    m_dot_reducer = vol_dot_press_fuel*rho_reducer_fuel 


    
    mass_press -= m_dot_reducer*dt
    
    rho_press = mass_press/vol_press
    
    p_current = psi('P','D',rho_press,'T',T_current,press)
    
    mass_press_array_jt[ii] = mass_press
    T_reducer_fuel_array_jt[ii] = T_reducer_fuel
    T_reducer_ox_array_jt[ii] = T_reducer_ox
    T_presstank_array_jt[ii] = T_current

    # rho_reducer_array[ii] = rho_reducer
    p_press_array_jt[ii] = p_current
    
    
    
    #############
    
    m_dot_press_fuel = vol_dot_press_fuel*rho_reducer_fuel
    m_dot_press_ox = vol_dot_press_ox*rho_reducer_ox
    rho_std = psi('D','P',1.013e5,'T',293,'N2')
    q_std_fuel = m_dot_press_fuel / rho_std    # std m^3/s
    q_std_ox = m_dot_press_ox / rho_std    
    q_std_fuel = q_std_fuel*1000*60 # std l/min
    q_std_ox = q_std_ox*1000*60
    
    q_std_fuel_list.append(q_std_fuel)
    q_std_ox_list.append(q_std_ox)
    
    cv_fuel.append(get_cv(q_std_fuel, p_current, T_current, 'N2'))
    cv_ox.append(get_cv(q_std_ox, p_current, T_current, 'N2'))

#%%
plt.figure(figsize=(12,8))

plt.subplot(221)
plt.plot(np.linspace(0, deltat, N), T_presstank_array_jt, label='bundle')
plt.plot(np.linspace(0, deltat, N), T_reducer_fuel_array_jt, label='fuel')
plt.plot(np.linspace(0, deltat, N), T_reducer_ox_array_jt, label='ox')
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.legend()
plt.grid()

plt.subplot(224)
plt.plot(np.linspace(0, deltat, N), cv_fuel, label='fuel')
plt.plot(np.linspace(0, deltat, N), cv_ox, label='ox')
plt.xlabel('Time [s]')
plt.ylabel('Cv value')
plt.legend()
plt.grid()

plt.subplot(222)
plt.plot(np.linspace(0, deltat, N), q_std_fuel_list, label='fuel')
plt.plot(np.linspace(0, deltat, N), q_std_ox_list, label='ox')
plt.xlabel('Time [s]')
plt.ylabel('Std l/min')
plt.legend()
plt.grid()

plt.subplot(223)
plt.plot(np.linspace(0, deltat, N), p_press_array_jt/1e5)#, label='fuel')
# plt.plot(np.linspace(0, deltat, N), q_std_ox_list, label='ox')
plt.xlabel('Time [s]')
plt.ylabel('Pressure [bar]')
# plt.legend()
plt.grid()