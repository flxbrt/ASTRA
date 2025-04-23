# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:09:13 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi


def get_kv(rho, deltap, Q):
    '''
    rho : TYPE
        kg/m^3.
    deltap : TYPE
        bar.
    Q : TYPE
        m^3/h.
    '''
    rho0 = psi('D','P',1e5,'T',293,'H2O') # always water!!!

    return Q*np.sqrt(rho/deltap/rho0)

def get_deltap(rho, kv, Q):
    
    rho0 = psi('D','P',1e5,'T',293,'H2O') # always water!!!

    return Q**2/kv**2*rho/rho0


p = 40e5
T = 90
fluid = 'O2'

rho = psi('D','P',p,'T',T,fluid)


m_dot = 0.4
Q_max = m_dot/rho*3600



# Q_std = 3.2

# Q_max = Q_std/1000/60*3600



Kv = 30/1000*60#1.8 #1e-3*60*30
# Kv = 1.6


# m_dot_min = 0.6*1/(1+0.9) #  ROF 0.9 at throttled conditions
# Q_min = m_dot_min/rho*3600


deltap_max = get_deltap(rho, Kv, Q_max)
# deltap_min = get_deltap(rho, 1.6, Q_min)




#%% gaseous



def get_kv_gaseous(T, p1, p2, Q):
    '''
    rhoN : TYPE
        kg/m^3.
    deltap : TYPE
        bar.
    Q : TYPE
        Nm^3/h.
    '''
    
    rhoN = psi('D','P',1e5,'T',293,'N2')
    
    if p2 > p1/2:
        return Q/514*np.sqrt(rhoN*T/(p1-p2)/p2)
    else:
        print('chocked')
        return Q/257/p1*np.sqrt(rhoN*T)
        
    


def get_deltap_gaseous():
    return


m_dot = 0.34                     # kg/s

rhoN = psi('D','P',1e5,'T',293,'N2')            # kg/m^3

Q_std = m_dot/rhoN                              # std m^3 / s
Q_std = Q_std*1000*60                           # std l / min

# Q_std = 50000                                    # std l / min

Q_max = Q_std/1000/60*3600 # [m^3/h]

T = 300
p1 = 300
p2 = 10
print(get_kv_gaseous(T, p1, p2, Q_max)/0.865, Q_max)