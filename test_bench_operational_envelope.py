# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 20:33:04 2024

@author: felix
"""



import cantera as ct
import numpy as np
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI as psi

def set_combustion():
    comb = ct.Solution('gri30_WARR.yaml')
    return comb

def obtain_combustion_properties(comb, p, fuel, fuel_cp, ROF, injec_t_fuel, injec_t_ox):
    def heat_capacity_ratio(comb):
        return comb.cp_mass / comb.cv_mass
    
    def get_delta_h(injec_t, ref_t, p, fluid):
        
        vap_t = psi('T','P',p,'Q',1,fluid)
        
        # ref_t: temperature with which cantera object is initialized
        # vap_t: vaporization temperature of fluid at given pressure
        # injec_t: propellant temperature when injected / tank propellant temperature
        # assumption: propellants are always liquid --> vaporisation temperature is always higher than injection temperature
        # 3 cases need to be distinguished
        # 1: injec_t < vap_t < ref_t
        if ref_t > vap_t:
            h_low = psi('H', 'P', p,'T', injec_t, fluid)
            h_high = psi('H', 'P', p,'T', ref_t, fluid)
            deltah = h_high - h_low
            # print(1)
        # 2: injec_t < ref_t < vap_t
        if ref_t < vap_t and ref_t > injec_t:
            deltah = psi('H', 'P', p,'T', ref_t, fluid) - psi('H', 'P', p,'T', injec_t, fluid) + psi('H','P',p,'Q',1,fluid) - psi('H','P',p,'Q',0,fluid)
            # print(2)
        # 3: ref_t < injec_t < vap_t
        if ref_t < injec_t:
            # negative sign for the first paranthesis, since deltah is subtracted and therefore negative frist term yields less enthalpy being subtracted from cantera object
            deltah =  - (psi('H', 'P', p,'T', injec_t, fluid) - psi('H', 'P', p,'T', ref_t, fluid)) + psi('H','P',p,'Q',1,fluid) - psi('H','P',p,'Q',0,fluid)
            # print(3)
        return deltah 
    
    ref_t = 273

    comb.Y = {fuel: 1, 'O2': ROF}
    comb.TP = ref_t, p
    comb.equilibrate('HP')
    
    delta_h_ox = get_delta_h(injec_t=injec_t_ox, ref_t=ref_t, p=p, fluid='O2')
    delta_h_fuel = get_delta_h(injec_t=injec_t_fuel, ref_t=ref_t, p=p, fluid=fuel_cp)
    
    deltah = ROF/(1+ROF)*delta_h_ox + 1/(1+ROF)*delta_h_fuel
    
    comb.HP = comb.h-deltah, p
    
    return comb.T, comb.mean_molecular_weight, heat_capacity_ratio(comb)

def calc_chocked_mass_flow(d_th, gamma, M, T, p_cc):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        m_dot = np.pi*d_th**2/4*p_cc/np.sqrt(T*R_ideal/M)*theta
        return m_dot
    
def calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M):
    def expansion_ratio(p_e, p_cc, gamma):
        nominator = ((gamma-1)/2)*(2/(gamma+1))**((gamma+1)/(gamma-1))
        denominator = (p_e/p_cc)**(2/gamma)*(1-(p_e/p_cc)**((gamma-1)/gamma))
        return np.sqrt(nominator/denominator)
    
    def pressure_ratio(epsilon, gamma):
        def f(p_ratio, epsilon, gamma):
            nominator = ((gamma-1)/2)*(2/(gamma+1))**((gamma+1)/(gamma-1))
            denominator = (p_ratio)**(2/gamma)*(1-(p_ratio)**((gamma-1)/gamma))
            return epsilon - np.sqrt(nominator/denominator)
        p_ratio = fsolve(f, 0.01, args=(epsilon, gamma))
        return p_ratio 
        
    
    
    
    p_ratio = pressure_ratio(epsilon, gamma)
    
    p_exit = p_cc*p_ratio
    
    A_exit = d_th**2/4*np.pi*epsilon
    
    R_ideal = 8314
    
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_exit/p_cc)**((gamma-1)/gamma)))

    p_amb = 1e5
    
    v_eff = v_ex + A_exit/m_dot*(p_exit-p_amb)
    
    return v_eff
    
comb = set_combustion()
fuel = 'C2H5OH'
fuel_cp = 'C2H6O'
p_cc = 30e5
ROF = 1.025
injec_t_fuel = 400
injec_t_ox = 100
T, M, gamma = obtain_combustion_properties(comb, p_cc, fuel, fuel_cp, ROF, injec_t_fuel, injec_t_ox)

d_th = 0.0303
epsilon = 2.92447450715475

nozzle_eff = 0.9
comb_eff = 0.9
eff = nozzle_eff*comb_eff

# eff = 0.912

#%% case SM1

#!!! you have to change efficiencies above!!!

# m_dot = calc_chocked_mass_flow(d_th, gamma, M, T, p_cc)

# v_eff = calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M)[0]

# v_eff = v_eff*eff # 90% comb efficiency

# thrust = m_dot * v_eff

# m_dot_ox = m_dot * ROF/(1+ROF)
# m_dot_fuel = m_dot * 1/(1+ROF)

#%% case SM2

#!!! you have to change efficiencies above!!!

# m_dot = calc_chocked_mass_flow(d_th, gamma, M, T, p_cc)

# v_eff = calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M)[0]

# v_eff = v_eff*eff # 90% comb efficiency

# thrust = m_dot * v_eff

# m_dot_scaled = m_dot*2.75e3/thrust # scaled towards 2.7kN

# m_dot_ox = m_dot_scaled * ROF/(1+ROF)
# m_dot_fuel = m_dot_scaled * 1/(1+ROF)

#%% case SW1 --> using values from before

# m_dot = calc_chocked_mass_flow(d_th, gamma, M, T, p_cc)

# v_eff = calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M)[0]

# v_eff = v_eff*eff # 90% comb efficiency

# thrust = m_dot * v_eff

# m_dot_scaled = m_dot*3e3/thrust # scaled towards 3kN

# m_dot_ox = m_dot_scaled * ROF/(1+ROF)
# m_dot_fuel = m_dot_scaled * 1/(1+ROF)

#%%  case TW2

#!!! you have to change LOx injection temperature above!!!

# m_dot = calc_chocked_mass_flow(d_th, gamma, M, T, p_cc)

# v_eff = calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M)[0]

# v_eff = v_eff*eff # 95% comb efficiency

# thrust = m_dot * v_eff

# m_dot_scaled = m_dot*2.75e3/thrust # scaled towards 2.7kN

# m_dot_ox = m_dot_scaled * ROF/(1+ROF)
# m_dot_fuel = m_dot_scaled * 1/(1+ROF)

#%% case RM1 & RW1

#!!! you have to change ROF above!!!

m_dot = calc_chocked_mass_flow(d_th, gamma, M, T, p_cc)

v_eff = calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M)[0]

v_eff = v_eff*eff # 90% comb efficiency

thrust = m_dot * v_eff

m_dot_scaled = m_dot*2.75e3/thrust # scaled towards 2.7kN

m_dot_ox = m_dot_scaled * ROF/(1+ROF)
m_dot_fuel = m_dot_scaled * 1/(1+ROF)

#%% case PW1

#!!! you have to change p_cc above!!!

# m_dot = calc_chocked_mass_flow(d_th, gamma, M, T, p_cc)

# v_eff = calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M)[0]

# v_eff = v_eff*eff # 90% comb efficiency

# thrust = m_dot * v_eff

# m_dot_scaled = m_dot*2.75e3/thrust # scaled towards 2.7kN

# m_dot_ox = m_dot_scaled * ROF/(1+ROF)
# m_dot_fuel = m_dot_scaled * 1/(1+ROF)