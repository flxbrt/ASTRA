# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 07:28:33 2026

@author: felix
"""



import yaml
import numpy as np
from ignitor import main
from volume_massflow_conversion import continuity
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as psi



def plot_over_rof(rof, T_ign, c_star, m_fu, m_ox, V_fu, V_ox,
                  d_th, d_fu, d_ox, d_q, mie, V_bd, V_fu_tank, V_ox_tank, v_fu, v_ox, T_max = None,
                  flb_lim=False):

    # flammability limits
    rho_h2 = psi('D', 'P', 1e5 ,'T', 288, 'H2')
    rho_air = psi('D', 'P', 1e5 ,'T', 288, 'Air')
    
    v_lfl_h2 = 0.04
    v_ufl_h2 = 0.75
    
    v_lfl_air = 1 - v_lfl_h2
    v_ufl_air = 1 - v_ufl_h2

    rof_lfl = v_lfl_air / v_lfl_h2 * rho_air / rho_h2
    rof_ufl = v_ufl_air / v_ufl_h2 * rho_air / rho_h2

    def add_fl_limits(ax):
        if flb_lim:
            ax.axvline(x=rof_lfl, label='LFL', color='red')
            ax.axvline(x=rof_ufl, label='UFL', color='purple')
            ax.legend()

    fig, axs = plt.subplots(4, 3, figsize=(18, 12))

    # 1) Combustion temperature
    axs[0, 0].plot(rof, T_ign)
    axs[0, 0].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[0, 0].set_ylabel("Combustion Temperature (K)")
    add_fl_limits(axs[0, 0])
    axs[0, 0].grid(True)

    # 2) c* / T_wall
    if T_max is not None:
        axs[1, 0].plot(rof, T_max)
        axs[1, 0].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
        axs[1, 0].set_ylabel("Wall Temperature at End of Burn (K)")
        add_fl_limits(axs[1, 0])
        axs[1, 0].grid(True)
    else:
        axs[1, 0].plot(rof, c_star)
        axs[1, 0].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
        axs[1, 0].set_ylabel("c* (m/s)")
        add_fl_limits(axs[1, 0])
        axs[1, 0].grid(True)

    # 3) Mass flows
    axs[0, 1].plot(rof, m_fu, label="Fuel")
    axs[0, 1].plot(rof, m_ox, label="Oxidizer")
    axs[0, 1].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[0, 1].set_ylabel("Mass Flow (g/s)")
    add_fl_limits(axs[0, 1])
    axs[0, 1].legend()
    axs[0, 1].grid(True)

    # 4) Volume flows
    axs[1, 1].plot(rof, V_fu, label="Fuel")
    axs[1, 1].plot(rof, V_ox, label="Oxidizer")
    axs[1, 1].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[1, 1].set_ylabel("Volume Flow (std l/min)")
    add_fl_limits(axs[1, 1])
    axs[1, 1].legend()
    axs[1, 1].grid(True)

    # 5) Throat diameter
    axs[0, 2].plot(rof, d_th)
    axs[0, 2].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[0, 2].set_ylabel("Throat Diameter (mm)")
    add_fl_limits(axs[0, 2])
    axs[0, 2].grid(True)

    # 6) Injector diameters
    axs[1, 2].plot(rof, d_fu, label="Fuel")
    axs[1, 2].plot(rof, d_ox, label="Oxidizer")
    axs[1, 2].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[1, 2].set_ylabel("Injector Diameter (mm)")
    add_fl_limits(axs[1, 2])
    axs[1, 2].legend()
    axs[1, 2].grid(True)

    # 7) Quenching distance
    axs[2, 0].plot(rof, d_q)
    axs[2, 0].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[2, 0].set_ylabel("Quenching Distance (mm)")
    add_fl_limits(axs[2, 0])
    axs[2, 0].grid(True)

    # 8) MIE
    axs[2, 1].plot(rof, mie)
    axs[2, 1].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[2, 1].set_ylabel("MIE (mJ)")
    add_fl_limits(axs[2, 1])
    axs[2, 1].grid(True)

    # 9) Breakdown voltage
    axs[2, 2].plot(rof, V_bd)
    axs[2, 2].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[2, 2].set_ylabel("Breakdown Voltage (kV)")
    add_fl_limits(axs[2, 2])
    axs[2, 2].grid(True)
    
    # 10) Storage volume for tanks
    axs[3, 0].plot(rof, V_fu_tank, label="Fuel")
    axs[3, 0].plot(rof, V_ox_tank, label="Oxidizer")
    axs[3, 0].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[3, 0].set_ylabel("Tank Volume (l)")
    add_fl_limits(axs[3, 0])
    axs[3, 0].legend()
    axs[3, 0].grid(True)
    
    # 11) Injection velocities
    axs[3, 1].plot(rof, v_fu, label="Fuel")
    axs[3, 1].plot(rof, v_ox, label="Oxidizer")
    axs[3, 1].set_xlabel("Oxidizer-to-Fuel Ratio (-)")
    axs[3, 1].set_ylabel("Injection Velocity (m/s)")
    add_fl_limits(axs[3, 1])
    axs[3, 1].legend()
    axs[3, 1].grid(True)

    plt.tight_layout()
    plt.show()



if __name__ == '__main__':

    rof_range = np.linspace(13, 102, 100) # 13 is the lowr limit and 102 is the upper limit for the applicability of the experimental interpolation
    
    # rof_range = [40]
    
    T_ign = list()
    
    c_star = list()
    
    m_fu = list()
    
    m_ox = list()
    
    d_th = list()
    
    V_ox = list()
    
    V_fu = list()
    
    d_fu = list()
    
    d_ox = list()
    
    d_q = list()
    
    V_bd = list()
    
    mie = list()
    
    T_max = list()
    
    V_fu_tank = list()
    
    V_ox_tank = list()
    
    v_fu = list()
    
    v_ox = list()
    
    # with open('chamber_warr.yaml', 'r') as f:
    with open('chamber.yaml', 'r') as f:
        chamber = yaml.load(f, Loader=yaml.SafeLoader)
    with open('ignitor.yaml', 'r') as f:
        ignitor_ipt = yaml.load(f, Loader=yaml.SafeLoader)
        
    fu = continuity(ignitor_ipt['fluid']['fu'])
    ox = continuity(ignitor_ipt['fluid']['ox'])
    unit = 'g/s'
    
    p = 1e5     # can be anything -> is anyway scaled to normal conditions
    T = 293
    
    for of in rof_range:
        
        print(f'Iter {np.where(rof_range==of)[0][0]+1} of {len(rof_range)}')
        
        ignitor_ipt['design']['of'] = of
        
        ignitor = main(chamber, ignitor_ipt)
        
        T_ign.append(ignitor['output']['T_ign'])
        c_star.append(ignitor['output']['c_star'])
        
        d_th.append(ignitor['output']['d_th']*1000)

        d_fu.append(ignitor['output']['d_fu']*1000)
        d_ox.append(ignitor['output']['d_ox']*1000)
        
        value = ignitor['output']['m_fu']*1000
        m_fu.append(value)
        fu.to_standard(value, unit)
        fu.qdot_mdot(p=p, T=T, from_to='m_to_q')
        fu.to_unit('std l/min')
        V_fu.append(fu.scaled_value)
        
        value = ignitor['output']['m_ox']*1000
        m_ox.append(value)
        ox.to_standard(value, unit)
        ox.qdot_mdot(p=p, T=T, from_to='m_to_q')
        ox.to_unit('std l/min')
        V_ox.append(ox.scaled_value)
        
        d_q.append(ignitor['output']['d_q']*1000)
        mie.append(ignitor['output']['MIE']*1000)
        V_bd.append(ignitor['output']['V_bd']/1000)
        
        T_max.append(ignitor['output']['T_max'])
        
        V_fu_tank.append(ignitor['output']['V_fu_tank']*1000)
        V_ox_tank.append(ignitor['output']['V_ox_tank']*1000)
        
        v_fu.append(ignitor['output']['v_fu'])
        v_ox.append(ignitor['output']['v_ox'])
        
    plot_over_rof(rof_range, T_ign, c_star, m_fu, m_ox, V_fu, V_ox, d_th, d_fu, d_ox, d_q, mie, V_bd, V_fu_tank, V_ox_tank, v_fu, v_ox, T_max)
