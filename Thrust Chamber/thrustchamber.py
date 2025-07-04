# -*- coding: utf-8 -*-
"""
Created on Sun May 25 11:07:32 2025

@author: felix
"""



import copy
import inspect
import numpy as np
import cantera as ct
from functools import partial
from scipy.optimize import root
from CoolProp.CoolProp import PropsSI as psi

class ThrustChamber:
    
    def __init__(self, fuel, ox, T_fuel, T_ox, pcc_nominal):
        self.fuel = fuel
        self.ox = ox
        self.pcc_nominal = pcc_nominal
        self.T_ox = T_ox
        self.T_fuel = T_fuel
        # self.rof = rof
        self.rho_fuel = psi('D', 'T', self.T_fuel, 'P', self.pcc_nominal, self.fuel) # assuming incompressible fluid
        self.rho_ox = psi('D', 'T', self.T_ox, 'P', self.pcc_nominal, self.ox) # assuming incompressible fluid
        
    def init_combustion(self):
        comb_base = ct.Solution('gri30_WARR.yaml')
        self.comb_base = ct.Quantity(comb_base, constant='HP')
        
    def set_combustion(self, p, rof):
        comb = copy.copy(self.comb_base)

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
        
        if self.fuel=='C2H6O':
            fuel_ct = 'C2H5OH'
        else:
            fuel_ct = self.fuel
        
        ref_t = 273
        comb.Y = {fuel_ct: 1, self.ox: rof}
        comb.TP = 273, p
        comb.equilibrate('HP')
        
        delta_h_ox = get_delta_h(injec_t=self.T_ox, ref_t=ref_t, p=p, fluid=self.ox)
        delta_h_fuel = get_delta_h(injec_t=self.T_fuel, ref_t=ref_t, p=p, fluid=self.fuel)
        deltah = rof/(1+rof)*delta_h_ox + 1/(1+rof)*delta_h_fuel
        comb.HP = comb.h-deltah, p
        self.comb = comb
        
        self.gamma = comb.cp_mass / comb.cv_mass
        self.Tcc = comb.T
        self.Mcc = comb.mean_molecular_weight
        R = 8314
        self.Rcc = R/self.Mcc
        self.Cpcc = comb.cp_mass
        self.lam = comb.thermal_conductivity
        self.Prcc = comb.viscosity*comb.cp_mass/comb.thermal_conductivity
        
    def get_c_star(self, eta=1):
        c_star = np.sqrt(self.Rcc*self.Tcc/self.gamma*((self.gamma+1)/2)**((self.gamma+1)/(self.gamma-1)))
        return c_star*eta
    
    def get_thrust_coefficient(self, m_dot, ve, A_th, eps_e, p_cc, p_e, beta=0, p_amb=1e5):
        eta = np.cos(np.deg2rad(beta))
        cF = m_dot*ve*eta/p_cc/A_th*eps_e*((p_e-p_amb)/p_cc)
        return cF
    
    def get_expansion_ratio(self, p_cc, p_e):
        gamma = self.gamma
        nominator = ((gamma-1)/2)*(2/(gamma+1))**((gamma+1)/(gamma-1))
        denominator = (p_e/p_cc)**(2/gamma)*(1-(p_e/p_cc)**((gamma-1)/gamma))
        return np.sqrt(nominator/denominator)

    def get_throat_diameter(self, m_dot, p_cc):
        gamma = self.gamma
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        # d_th = np.sqrt(m_dot/p_cc/np.pi*4/np.sqrt(M/(R_ideal*T))*theta)
        d_th = np.sqrt(4/np.pi*m_dot/p_cc*np.sqrt(self.Tcc*self.Rcc)/theta)
        return d_th
    
    def get_chocked_mass_flow(self, d_th, p_cc):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        gamma = self.gamma
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        m_dot = np.pi*d_th**2/4*p_cc/np.sqrt(self.Tcc*R_ideal/self.Mcc)*theta
        return m_dot
    
    def get_isentropic_exhaust_velocity(self, p_cc, p_ex):
        R_ideal = 8314
        gamma = self.gamma
        # comb = set_combustion()
        v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(self.Tcc)/(self.Mcc)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))
        return v_ex
    
    def get_effective_exhaust_velocity(self, m_dot, pcc, pe, d_th, eps_e, eta_nozzle, pamb=1e5):
        Ae = d_th**2/4*np.pi*eps_e
        ve = self.get_isentropic_exhaust_velocity(pcc, pe)
        v_eff = ve*eta_nozzle + Ae/m_dot*(pe-pamb)
        return v_eff
        
    def get_thrust(self, m_dot, ve, Ae, pe, pamb=1e5):
        return m_dot*ve + Ae*(pe - pamb) 

    def get_massflow_from_thrust(self, F, Ae, ve, pe, pamb=1e5):
        return (F - Ae*(pe - pamb)) / ve
        
    def solve_for_thrust_constraint(self, F_max, F_min, pe_min, rof_min, rof_max, eta_cc=1, beta=0, ret_ext=False):
        '''
        solves for the throat diameter and combustion chamber pressure 
        corresponding to the minimum thrust given the values of maximum 
        and minimum thrust
        high thrust load point determines throat diameter
        low thrust load point determines expansion ratio
        '''
        eta_nozzle = np.cos(np.deg2rad(beta))
        pcc_max = self.pcc_nominal
        
        ### ---------------------- initial guess values ----------------------- 
        pcc_min = pcc_max*F_min/F_max 
        pe_max = pe_min*pcc_max/pcc_min # pcc/pe is constant for all operating points; only holds true for isentropic flow and gamma being constant

        self.init_combustion()
        self.set_combustion(p=pcc_max, rof=rof_max)
        ve_max = self.get_isentropic_exhaust_velocity(pcc_max, pe_max)*eta_nozzle
        mdot_max = F_max / ve_max / eta_cc # rough approximate since ce not yet available
        d_th = self.get_throat_diameter(m_dot=mdot_max, p_cc=pcc_max)
        
        eps_e = self.get_expansion_ratio(p_cc=pcc_min, p_e=pe_min)
        d_e = np.sqrt(eps_e)*d_th
        A_e = d_e**2/4*np.pi
        
        F_max_comp = self.get_thrust(mdot_max, ve_max, A_e, pe_max)
        mdot_min = self.get_chocked_mass_flow(d_th=d_th, p_cc=pcc_min)/eta_cc
        ve_min = self.get_isentropic_exhaust_velocity(p_cc=pcc_min, p_ex=pe_min)*eta_nozzle
        F_min_comp = self.get_thrust(mdot_min, ve_min, A_e, pe_min)
        
        ### ------------------------- main iteration --------------------------
        err = 1
        while err > 0.0001:
            
            # compute min thrust
            self.set_combustion(p=pcc_min, rof=rof_min)
            ve_min = self.get_isentropic_exhaust_velocity(p_cc=pcc_min, p_ex=pe_min)
            ve_min_eta = ve_min * eta_nozzle# part of exhaust velcity that contributes to the thrust
            mdot_min = self.get_chocked_mass_flow(d_th=d_th, p_cc=pcc_min) # mass flow that would result from pcc_min if there was no combustion efficiency
            mdot_min_eta = mdot_min/eta_cc # considers extra mass flow required to set min combustion chamber pressure given combustion efficiency
            F_min_comp = self.get_thrust(mdot_min_eta, ve_min_eta, A_e, pe_min)
            # correct pressure to converge on min thrust required
            deltap = (F_min-F_min_comp)*1e5/100*0.5 # adaptive step size --> 1bar causes around 100N change --> step size is scaled such that 100N offset causes 0.5bar adaption and therefore slow but robust convergence
            pcc_min += deltap
            # compute expansion ratio
            eps_e = self.get_expansion_ratio(p_cc=pcc_min, p_e=pe_min)
            d_e = np.sqrt(eps_e)*d_th
            A_e = d_e**2/4*np.pi
            
            ce_min = self.get_effective_exhaust_velocity(mdot_min_eta, pcc_min, pe_min, d_th, eps_e, eta_nozzle)
            
            # compute max thrust
            pe_max = pe_min*pcc_max/pcc_min # pcc/pe is constant for all operating points; only holds true for isentropic flow and gamma being constant
            self.set_combustion(p=pcc_max, rof=rof_max)
            ve_max = self.get_isentropic_exhaust_velocity(pcc_max, pe_max)
            ve_max_eta = ve_max * eta_nozzle
            mdot_max_eta = self.get_massflow_from_thrust(F_max, A_e, ve_max_eta, pe_max) # this is the mass flow required to obtain the computed thrust INCLUDING combustion efficiency considerations
            mdot_max = mdot_max_eta * eta_cc # if there was a perfect combustion efficiency, only this mass flow was required
            F_max_comp = self.get_thrust(mdot_max_eta, ve_max_eta, A_e, pe_max)
            # compute throat diameter
            d_th = self.get_throat_diameter(m_dot=mdot_max, p_cc=pcc_max) # use mass flow for perfect combustion efficiency
            
            ce_max = self.get_effective_exhaust_velocity(mdot_max_eta, pcc_max, pe_max, d_th, eps_e, eta_nozzle)
            
            err = abs(F_min_comp - F_min)/F_min

        if ret_ext:
            return d_th, pcc_min, ce_min, ce_max, ve_min, mdot_min, mdot_min_eta, F_min_comp, pe_max, ve_max, mdot_max, mdot_max_eta, F_max_comp, eps_e, eta_nozzle
        else:
            return d_th, pcc_min, ce_min, ce_max
        
            
fuel = 'C2H6O'
ox = 'O2'
T_fuel = 400
T_ox = 110
pcc_nominal = 20e5
F_nom = 2834#2500
F_min = 1154#0.4*F_nom
# rof = 1
pe_min = 0.6e5
rof_min = 0.92
rof_max = 1.17
tc = ThrustChamber(fuel, ox, T_fuel, T_ox, pcc_nominal)#, rof)

d_th, pcc_min, ce_min, ce_max, ve_min, mdot_min, mdot_min_eta, F_min_comp, pe_max, ve_max, mdot_max, mdot_max_eta, F_max_comp, eps_e, eta_nozzle = tc.solve_for_thrust_constraint(
    F_max=F_nom, F_min=F_min, pe_min=pe_min, rof_min=rof_min, rof_max=rof_max, eta_cc=0.95, beta=15, ret_ext=True)

# print(f'Min effective Isp: {(ce_max+ce_min)/2:.2f}[m/s]')

# # isp im mittel ausrechnen
# # tabelle für matej machen
# eta = 1
# tc.set_combustion(p=pcc_min, rof=rof_min)
# c_star_min = tc.get_c_star(eta)

# tc.set_combustion(p=pcc_nominal, rof=rof_max)
# c_star_max = tc.get_c_star(eta)