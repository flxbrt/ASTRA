# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 10:54:36 2026

@author: felix
"""



import yaml
import numpy as np
import cantera as ct
from CoolProp.CoolProp import PropsSI as psi
from scipy.optimize import fsolve
from scipy.interpolate import interp1d




def heat_of_combustion(fuel, ox, of):
    
    if fuel=='H2':
        Hu = 120e6 # MJ/kg
        of_s = 8   # stoichiometric of for reaction with O2
    elif fuel=='C2H5OH' or fuel=='C2H6O':
        Hu = 26.8e6   # MJ/kg
        of_s = 2.09 # stoichiometric of for reaction with O2

    # correct for oxidizer of choice -> only relevant to ignitor and hence H2 as fuel
    if fuel == 'H2' and ox == 'Air':
        fac = 8/34.5 # stoichiometric ratio for H2-air is 34.5
    else:
        fac = 1
        
    of_s /= fac

    # compute heat release given rof
    if of >= of_s: # all fuel combusted
        return Hu
    else: # only part of the fuel is combusted
        return Hu*of/of_s
    


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
        # print(deltah)
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



def igniter_power(chamber, ignitor):
    
    key = ignitor['models']['key']    
    
    margin = ignitor['models']['P_s']        # margin for ignitor power added on top
    
    m_mc = chamber['m_t'] 
    
    ##### option 1: given thermal power released by main chamber (Max Finkl)
    if key == 'to':
                  # kg/s
        Hu = heat_of_combustion(chamber['fu_cp'], chamber['ox'], chamber['of'])
        P_mc_to = m_mc*Hu                  # W thermal heat release of main chamber
        
        P_ign_to = ignitor['models']['P_to']*P_mc_to    # W thermal input required by igniter
        
        P_ign_to *= (1 + margin)
        
        ignitor['output']['P_ign'] = P_ign_to
    
    ##### option 2: given temperature increase to autoignition temperature (Francesco Palmese)
    elif key == 'ai':
        T_mc_ai = chamber['T_ai']
    
        of_mc = chamber['of']
        m_mc_fu = 1/(1+of_mc)*m_mc
        m_mc_ox = of_mc/(1+of_mc)*m_mc
        
        # print(m_mc_ox/m_mc_fu) # senity check -> correct
        
        deltah_fu = get_delta_h(injec_t=chamber['T_fu'], ref_t=T_mc_ai, p=chamber['p_c'], fluid=chamber['fu_cp'])
        deltah_ox = get_delta_h(injec_t=chamber['T_ox'], ref_t=T_mc_ai, p=chamber['p_c'], fluid=chamber['ox'])
        
        P_mc_ai = deltah_fu * m_mc_fu + deltah_ox * m_mc_ox # W thermal input required to heat propellants to auto ignition temperature
    
        P_ign_ai = ignitor['models']['P_ai']*P_mc_ai    # W thermal input required by igniter
    
        P_ign_ai *= (1 + margin)
    
        ignitor['output']['P_ign'] = P_ign_ai
    
    else:
        raise ValueError(f'{key=} invalid')
    
    return ignitor
    

def ignitor_combustion(chamber, ignitor):
    
    of_ign = ignitor['design']['of']
    
    fu = ignitor['fluid']['fu']
    ox = ignitor['fluid']['ox']
    
    p_ign = chamber['p_c']*ignitor['design']['Pi_p']
    
    ignitor['output']['p_ign'] = p_ign
    
    
    if ox == 'Air':
        
        # Air: 23.2% O2, 76.8% N2 by mass (= ~21/79 by mole)
        air_o2_massfrac = 0.21
        air_n2_massfrac = 0.79
        ox_dict = {
            'O2': of_ign * air_o2_massfrac,
            'N2': of_ign * air_n2_massfrac
        }
        
    else:
        ox_dict = {ox: of_ign}
    
    composition = {fu: 1.0, **ox_dict}
    
    comb = ct.Solution('h2o2.yaml')              # obtain enthalpy at equilibrium conditions
    
    comb.Y = composition
    
    T_fu = ignitor['fluid']['T_fu']
    T_ox = ignitor['fluid']['T_ox']
    T_mean = (T_fu + T_ox*of_ign)/(1+of_ign)        # mass averaged temperature
    
    comb.TP = T_mean, p_ign
    comb.equilibrate('HP') 
    
    T_ign = comb.T
    M_ign = comb.mean_molecular_weight 
    gamma = comb.cp_mass / comb.cv_mass  

    c_star_id = np.sqrt(8314/M_ign*T_ign/gamma*((gamma+1)/2)**((gamma+1)/(gamma-1)))
    c_star_eff = c_star_id * ignitor['models']['eta']  # combusiton efficiency
    
    ignitor['output']['T_ign'] = T_ign
    ignitor['output']['M_ign'] = M_ign
    ignitor['output']['gamma'] = gamma
    ignitor['output']['c_star_id'] = c_star_id
    ignitor['output']['c_star'] = c_star_eff            # effective c star given combustion efficiency

    return ignitor, comb


def igniter_mass_flow(chamber, ignitor, comb):
    # option 1 and 2 refer to igniter power options
    
    of_ign = ignitor['design']['of']
    
    P_ign= ignitor['output']['P_ign']
    
    fu = ignitor['fluid']['fu']
    ox = ignitor['fluid']['ox']
    
    key = ignitor['models']['key']
    
    p_ign = ignitor['output']['p_ign']
    
    if key == 'to': ##### option 1 thermal power
        Hu = heat_of_combustion(fu, ox, of_ign)
        
        m_ign_fu = P_ign/Hu
        m_ign_ox = of_ign * m_ign_fu
    
    ##### option 2 auto ignition
    elif key == 'ai':
        if ox == 'Air':
            
            # Air: 23.2% O2, 76.8% N2 by mass (= ~21/79 by mole)
            air_o2_massfrac = 0.21
            air_n2_massfrac = 0.79
            ox_dict = {
                'O2': of_ign * air_o2_massfrac,
                'N2': of_ign * air_n2_massfrac
            }
            
        else:
            ox_dict = {ox: of_ign}
        
        composition = {fu: 1.0, **ox_dict}
        
        comb_ai = ct.Solution('h2o2.yaml')              # obtain enthalpy at autoignition conditions
        
        # mass averaged temperature
        
        T_ai = ignitor['models']['T_ai']
        
        comb_ai.Y = composition
        
        comb_ai.TP = T_ai, p_ign
        
        comb_ai.equilibrate('TP')
        
        h_ai = comb_ai.h
        
        comb_eq = comb
        
        # comb_eq.Y = composition
        
        # T_fu = ignitor['fluid']['T_fu']
        # T_ox = ignitor['fluid']['T_ox']
        # T_mean = (T_fu + T_ox*of_ign)/(1+of_ign)        # mass averaged temperature
        
        # comb_eq.TP = T_mean, p_ign
        # comb_eq.equilibrate('HP') # not required since enthalpy does not change
        
        h_eq = comb_eq.h
        
        deltah = h_eq - h_ai
        
        m_ign = P_ign / deltah
        
        m_ign_fu = m_ign * 1/(1+of_ign)
        m_ign_ox = m_ign * of_ign/(1+of_ign)
    
    # ideal mass flows
    ignitor['output']['m_fu_id'] = m_ign_fu#/eta
    ignitor['output']['m_ox_id'] = m_ign_ox#/eta
    ignitor['output']['m_tot_id'] = (m_ign_fu + m_ign_ox)#/eta
    
    # required mass flows to provide the required ignitor power at reduced combustion efficiency
    eta = ignitor['models']['eta']          # combusiton efficiency

    ignitor['output']['m_fu'] = m_ign_fu/eta
    ignitor['output']['m_ox'] = m_ign_ox/eta
    ignitor['output']['m_tot'] = (m_ign_fu + m_ign_ox)/eta
    
    return ignitor

def compute_geometry(ignitor):
    
    # def get_throat_diameter(gamma, m_dot, p, T, M):
        
    #     R = 8314/M
    #     # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
    #     theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
    #     # d_th = np.sqrt(m_dot/p_cc/np.pi*4/np.sqrt(M/(R_ideal*T))*theta)
    #     d_th = np.sqrt(4/np.pi*m_dot/p*np.sqrt(T*R)/theta)
    #     return d_th
    
    # p_ign = ignitor['output']['p_ign'] 
    # T_ign = ignitor['output']['T_ign'] 
    # M_ign = ignitor['output']['M_ign'] 
    # gamma = ignitor['output']['gamma']
    
    # dth = get_throat_diameter(gamma, m_tot, p_ign, T_ign, M_ign)
    
    def get_throat_diameter(ignitor, key=''):
        
        m_tot = ignitor['output'][f'm_tot{key}']
        
        c_star = ignitor['output']['c_star']    # efficiency already considered in m_tot, so effective c_star has to be called
        
        p_ign = ignitor['output']['p_ign']
        
        Ath = c_star*m_tot/p_ign
        dth = np.sqrt(4*Ath/np.pi)
        return dth
        
    dth_id = get_throat_diameter(ignitor, '_id')
    dth_eff = get_throat_diameter(ignitor)
    
    ignitor['output']['d_th_id'] = dth_id
    ignitor['output']['d_th'] = dth_eff

    return ignitor


def injector_cross_section(ignitor):
    
    # assuming isentropic compressible flow -> entire dynamic pressure is lost -> cd not considered
    
    cd = ignitor['models']['cd']
    
    T_fu = ignitor['fluid']['T_fu']
    T_ox = ignitor['fluid']['T_ox']
    
    fu = ignitor['fluid']['fu']
    ox = ignitor['fluid']['ox']
    
    dp_rel = ignitor['design']['Pi_deltap']
    p_ign = ignitor['output']['p_ign']          # downstream injector / static pressure
    p_inj = p_ign * (1 + dp_rel)                # upstream injector / total pressure
    
    ignitor['output']['p_inj'] = p_inj

    gamma_fu = psi('Cpmass', 'P', p_inj ,'T', T_fu, fu) / psi('Cvmass', 'P', p_inj ,'T', T_fu, fu)
    gamma_ox = psi('Cpmass', 'P', p_inj ,'T', T_ox, ox) / psi('Cvmass', 'P', p_inj ,'T', T_ox, ox)
    
    f = lambda Ma, gamma, p_ign, p_inj: p_inj/p_ign - (1+((gamma-1)/2)*Ma**2)**(gamma/(gamma-1))    # p_ign is static pressure / downstream of injector, p_inj is total pressure / upstream of injector
    
    Ma_fu = fsolve(f, x0=0.5, args=(gamma_fu, p_ign, p_inj))
    Ma_ox = fsolve(f, x0=0.5, args=(gamma_ox, p_ign, p_inj))
    
    T_fu_st = T_fu*(1 + ((gamma_fu-1)/2)*Ma_fu**2)**(-1)
    T_ox_st = T_ox*(1 + ((gamma_ox-1)/2)*Ma_ox**2)**(-1)
    
    # print(T_fu_st, T_ox_st)
    
    a_fu =  psi('speed_of_sound', 'P', p_ign ,'T', T_fu_st, fu)        # static mach number
    a_ox =  psi('speed_of_sound', 'P', p_ign ,'T', T_ox_st, ox)        # static mach number
    
    v_fu = a_fu*Ma_fu
    v_ox = a_ox*Ma_ox
    
    ignitor['output']['v_fu'] = v_fu
    ignitor['output']['v_ox'] = v_ox
    
    # previous version: assuming incompressible flow
    
    # dp = p_ign * dp_rel
    
    # rho_fu = psi('D', 'P', p_ign ,'T', T_fu, fu)
    # rho_ox = psi('D', 'P', p_ign ,'T', T_ox, ox)
    
    # v_fu = np.sqrt(2*dp/rho_fu/cd)
    # v_ox = np.sqrt(2*dp/rho_ox/cd)
    
    # a_fu =  psi('speed_of_sound', 'P', p_inj ,'T', T_fu, fu)
    # a_ox =  psi('speed_of_sound', 'P', p_inj ,'T', T_fu, ox)
    
    # if v_fu/a_fu > 0.3:
    #     print(f'Ma<0.3 assumption does not hold for {fu}')
    # if v_ox/a_ox > 0.3:
    #     print(f'Ma<0.3 assumption does not hold for {ox}')

    m_fu = ignitor['output']['m_fu']
    m_ox = ignitor['output']['m_ox']
    
    rho_fu = psi('D', 'P', p_ign ,'T', T_fu_st, fu)
    rho_ox = psi('D', 'P', p_ign ,'T', T_ox_st, ox)
    
    A_fu = m_fu/rho_fu/v_fu/cd
    A_ox = m_ox/rho_ox/v_ox/cd
    
    d_fu = np.sqrt(4*A_fu/np.pi)
    d_ox = np.sqrt(4*A_ox/np.pi)
    
    # print(d_fu)
    # print(d_ox)
    
    ignitor['output']['d_fu'] = d_fu
    ignitor['output']['d_ox'] = d_ox
    
    return ignitor

def spark_ignition(ignitor):
    
    def rof_to_conc(rof):
        ''' conversion from rof to h2 concentration 
        in air given partial pressure 
        which is identical to volumentric concentration for ideal gases '''
        
        rho_fu = psi('D', 'P', 1e5,'T', 288, 'H2')
        rho_ox = psi('D', 'P', 1e5,'T', 288, 'Air')
        
        return (rho_ox/rho_fu)/(rof + rho_ox/rho_fu)*100
    
    def minimum_ignition_energy(ignitor, verbose=False):
        
        ''' conservative estimate
        holds for any gap distance until 3mm
        '''
        of_ign = ignitor['design']['of']
        
        data = np.load(ignitor['data']['mie'])
        if verbose:
            print(f'Using {ignitor["data"]["mie"]} for MIE')
        
        h2 = data[0,:]
        mie = data[1,:]
        
        itp = interp1d(h2, mie, kind='linear')
        conc = rof_to_conc(of_ign)
        
        mie = itp(conc)
    
        return mie # in J
    
    def quenching_distance(ignitor, verbose=False):
        
        of_ign = ignitor['design']['of']
        
        data = np.load(ignitor['data']['dq'])
        if verbose:
            print(f'Using {ignitor["data"]["dq"]} for quenching distance')
        
        h2 = data[0,:]
        dq = data[1,:]
        
        itp = interp1d(h2, dq, kind='linear')
        conc = rof_to_conc(of_ign)
        
        dq = itp(conc)
        
        return dq # in m
    
    def breakdown_voltage(ignitor, verbose=False):
        
        of_ign = ignitor['design']['of']
        
        data = np.load(ignitor['data']['Vbd'])
        if verbose:
            print(f'Using {ignitor["data"]["Vbd"]} for breakdown voltage')
        
        h2 = data[0,:]
        Vbd = data[1,:]
        
        itp = interp1d(h2, Vbd, kind='linear')
        conc = rof_to_conc(of_ign)
        
        Vbd = itp(conc)
        
        return Vbd # in V

    verbose = ignitor['data']['verbose']

    mie = minimum_ignition_energy(ignitor, verbose)
    dq = quenching_distance(ignitor, verbose)
    Vbd = breakdown_voltage(ignitor, verbose)
    
    ignitor['output']['MIE'] = mie
    ignitor['output']['d_q'] = dq
    ignitor['output']['V_bd'] = Vbd
    
    ignitor['data']['verbose'] = 0
    
    return ignitor

def wall_temperature(ignitor, comb):
    
    def rk4(f, y, h, t, *args):
        # runge kutte 4th order explicit
        tk_05 = t + 0.5*h
        yk_025 = y + 0.5 * h * f(t, y, *args)
        yk_05 = y + 0.5 * h * f(tk_05, yk_025, *args)
        yk_075 = y + h * f(tk_05, yk_05, *args)
        
        return y + h/6 * (f(t, y, *args) + 2 * f(tk_05, yk_025, *args) + 2 * f(tk_05, yk_05, *args) + f(t+h, yk_075, *args))

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
    
    
    def get_hg(comb, d, T_w, m_dot):
        ''' Bartz equation implemented '''
        A = np.pi*d**2/4
        p0 = comb.P
        T0 = comb.T
        gamma = comb.cp / comb.cv
        a = np.sqrt(gamma * ct.gas_constant * comb.T / comb.mean_molecular_weight)
        mu = comb.viscosity 
        lam = comb.thermal_conductivity
        Pr = (mu * cp) / lam
        c_star = p0*d**2/4*np.pi/m_dot
        vel = m_dot / A / comb.density
        # print(vel)
        Ma = vel / a
        # print(Ma)
        C = 0.026
        omega = 0.6
        sigma = 1 / (((((1 / 2) * (T_w / T0) * (1 + (((gamma - 1) / 2) * Ma**2))) + (1 / 2))**(0.8 - (omega / 5))) * ((1 + (((gamma - 1) / 2) * Ma**2))**(omega / 5)))
        hg = (C / (d**0.2)) * ((mu**0.2 * cp)/(Pr**0.6)) * ((p0 / c_star)**0.8) * sigma # * fac
        
        return hg
    
    # Initial Conditions and Integration Parameters
    T_init = 300                # Initial Temperature
    t_start = 0                 # Simulation start time
    t_end = ignitor['design']['dur']  # Simulation end time
    h = 5e-5                    # time step size
    
    # Material properties
    rho = ignitor['design']['rho']
    cp = ignitor['design']['cp']
    lam = ignitor['design']['lam']
    
    # Ignitor properties
    m_tot = ignitor['output']['m_tot']
    d_th = ignitor['output']['d_th']

    # Geometry    
    r_inner = ignitor['output']['d_th']/2
    s = ignitor['design']['s']
    r_outer = r_inner + s
    l_unit = 1          # unit length of cylinder --> can be set to arbitrary value --> only required as senity checks to match physical units
    dr = 1e-4           # spatial discretisation
    radial_c = np.arange(r_inner, r_outer + 1e-6, dr) # radius at respective node centre
    radial_e = np.arange(r_inner-dr/2, r_outer+dr/2 + 1e-6, dr) # radius at node edge / interface between two nodes

    A = 2*np.pi*radial_e*l_unit     # cross section perpendicular to heat flux at each node edge
    V = (radial_e[1:]**2-radial_e[:-1]**2)*np.pi*l_unit # volume of each node / element
    m = rho*V                       # mass of each element

    # Numerical stability condition    
    condition = lam * h / (rho * cp * dr**2)
    if condition > 0.5:
        raise ValueError('Integration condition not satisfied')        
    
    time = np.arange(t_start, t_end, h)

    # Solution arrays
    temperature_array = np.ones((len(time)+1, len(radial_c)))*T_init
    internal_heat_array = np.ones((len(time)+1, len(radial_c)))*T_init*cp*m

    # Main Loop
    y = np.array([internal_heat_array[0,:], temperature_array[0,:]])

    for t_index, t in enumerate(time):
        T_w = y[1,0]
        
        hg = get_hg(comb, d_th, T_w, m_tot)
        
        Q_i = hg*(comb.T-T_w)*A[0]                    # inner BC heat flux into the wall; simplified assumption: stagnation instead of adiabatic wall temperature
        Q_e = 0                                     # outer BC heat flux exiting the wall
        
        args = [dr, cp, lam, m, Q_i, Q_e, A]
        y = rk4(time_derivative, y, h, t, *args)    # time integration of all radial nodes at once
        
        internal_heat_array[t_index+1, :] = y[0,:]
        temperature_array[t_index+1, :] = y[1,:]
        
    ignitor['output']['T_max'] = temperature_array[-1,0]
    
    return ignitor

def tank_volume(ignitor):
    
    p_stor = ignitor['design']['p_stor']
    T_fu = ignitor['fluid']['T_fu']
    T_ox = ignitor['fluid']['T_ox']
    fu = ignitor['fluid']['fu']
    ox = ignitor['fluid']['ox']

    dur = ignitor['design']['dur']                  # burn duration
    
    m_fu = ignitor['output']['m_fu']
    m_ox = ignitor['output']['m_ox']

    rho_fu = psi('D', 'P', p_stor,'T', T_fu, fu)    # storage density
    rho_ox = psi('D', 'P', p_stor,'T', T_ox, ox)
    
    V_fu_tank = m_fu*dur / rho_fu                   # ideal storage volume -> optimistic
    V_fu_ox = m_ox*dur / rho_ox
    
    ignitor['output']['V_fu_tank'] = V_fu_tank
    ignitor['output']['V_ox_tank'] = V_fu_ox

    return ignitor

def main(chamber, ignitor, eval_T=True):
    
    ignitor['output'] = {}
    
    ignitor = igniter_power(chamber, ignitor)
    
    ignitor, comb = ignitor_combustion(chamber, ignitor)
    
    ignitor = igniter_mass_flow(chamber, ignitor, comb)
    
    ignitor = compute_geometry(ignitor)

    ignitor = injector_cross_section(ignitor)
    
    ignitor = spark_ignition(ignitor)
    
    if eval_T:
    
        ignitor = wall_temperature(ignitor, comb)
        
    ignitor = tank_volume(ignitor)
    
    return ignitor

#%%
if __name__ == '__main__':
    
    with open('chamber.yaml', 'r') as f:
        chamber = yaml.load(f, Loader=yaml.SafeLoader)
    with open('ignitor.yaml', 'r') as f:
        ignitor = yaml.load(f, Loader=yaml.SafeLoader)
    
    ignitor = main(chamber, ignitor)