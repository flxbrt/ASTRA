# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 10:33:19 2025

@author: felix
"""



import numpy as np



class screw():
    
    def __init__(self):
        # Tabellenbuch Metall p. 214
        self.core_diameters = {
        "M4": 3.14,
        "M5": 4.02,
        "M6": 4.77,
        "M8": 6.47,
        "M10": 8.16,
        "M12": 9.85,
        "M14": 11.55
        }
        
        self.cross_section = {
        "M4": 8.78,
        "M5": 14.2,
        "M6": 20.1,
        "M8": 36.6,
        "M10": 58.0,
        "M12": 84.3,
        "M14": 115.47
        }

        self.yield_strength = {
        "5.8": 400,
        "6.8": 480,
        "8.8": 640,
        "9.8": 720,
        "10.9": 900,
        "12.9": 1080
        }
        
        self.flanken_d = {
        "M4": 3.55,
        "M5": 4.48,
        "M6": 5.35,
        "M8": 7.19,
        "M10": 9.03,
        "M12": 10.86,
        "M14": 12.7
        }
        
        self.steigung = {
        "M4": 0.7,
        "M5": 0.8,
        "M6": 1,
        "M8": 1,
        "M10": 1.5,
        "M12": 1.75,
        "M14": 2
        }
        
        self.d_W = {
        "M4": 5.9,
        "M5": 6.9,
        "M6": 8.9,
        "M8": 11.6,
        "M10": 14.6
        }
        
    def comp_tension(self, F, size):
        # d = self.core_diameters[size]*1e-3
        # A = d**2/4*np.pi
        A = self.cross_section[size]*1e-6
        return F/A
    
    def find_strength_class(self, sigma):
        for key, strength in sorted(self.yield_strength.items(), key=lambda x: x[1]):
            if sigma <= strength*1e6:
                return key
        return "Keine Standardklasse ausreichend"
    
    def mounting_torque(self, F_pretension, size, mu_g, mu_k):
        # senity checked with example from tabellenbuch --> yields the 127Nm
        d2 = self.flanken_d[size]*1e-3
        # print(d2)
        P = self.steigung[size]*1e-3
        phi = np.arctan(P/d2/np.pi)
        # print(np.rad2deg(phi))
        rho = np.arctan(mu_g/np.cos(np.deg2rad(60/2))) # metricshe gewinde haben 60° --> s. Zeichnung tabellenbuch S. 214
        # print(np.rad2deg(rho))
        
        dW = self.d_W[size]*1e-3 # Kopfdurchmesser
        Dki = int(size.split('M')[1])*1e-3 # S. 77 VDI 2230 Blatt 1
        dk = (dW+Dki)/2
        # print(dk)
        # alternativ Tabellenbuch S. 231
        # dk = 2*0.65*DKi
        
        # senity check from tabellenbuch metall p. 231
        # d2 = 10.86e-3
        # phi = np.deg2rad(2.9)
        # rho = np.deg2rad(10.5)
        # dk = 15.6e-3
        
        
        M = F_pretension*(d2/2*np.tan(phi+rho) + mu_k*dk/2)
        return M

#%%

################################### Ignitor ##################################
# screw itself has to carry the tension but also component in which screw is mounted has to carry tension

N_screws = 4

safety_fluid = 2
safety_screw = 2 # lower since you anyway round up; p.230 also defines to be >=1.5
mounting = 2

## Required pre tension --> only axially w.r.t. screw axis
# Force given sealing
r_groove_inner = 29*1e-3/2
width_groove = 2e-3
r_groove_outer = r_groove_inner + width_groove 

chord_thickness = 1.5e-3
groove_height = 1.1e-3

compression = (chord_thickness-groove_height)/chord_thickness

# approach 1: via e module
A_seal = np.pi*(r_groove_outer**2 - r_groove_inner**2)

# A_seal = np.pi*2*r_groove_inner*width_groove
el_mod = 17e6 # https://www.dt-bremen.de/UserFiles/file/pdf/Wissen/Wi-S02-06.pdf

sigma = compression*el_mod

F_seal = sigma*A_seal

# # approach 2; via parker handbook p. 114 --> Verformungskraft
# if compression < 0.1:
#     F_per_circum = 70*100 # Kraft pro Umfang [N/m]
# elif compression < 0.2:
#     F_per_circum = 190*100 # Kraft pro Umfang [N/m]
# elif compression < 0.3:
#     F_per_circum = 450*100 # Kraft pro Umfang [N/m]
# else:
#     raise ValueError('Compression exceeds limit')
    
# l_circum = r_groove_outer*2*np.pi

# F_seal = l_circum*F_per_circum

# Force given internal pressure
A_pressurised = r_groove_outer**2*np.pi
p_internal = 20e5

F_fluid = A_pressurised*p_internal

# overwriting the previous computation
# F_seal = 7000

# thermal expansion
# deltaT = 600
# alpha = 18e-6 # https://www.zollern.com/fileadmin/data/ZGM_Schmiede/files/de/datenblaetter/ZOLLERN-20200_Datenblatt_Schmiede_CCZr.pdf
# emod = 130e9
# A = 0.05**2/4*np.pi

# F_therm = A*emod*alpha*deltaT

# pure loading
F_load = F_fluid + F_seal

print(f'{F_seal=}\n{F_fluid=}\n{F_load=}')

screw_lox = screw()

# required pre tension with safety factor for fluid pressure
F_pretension = F_fluid*safety_fluid + F_seal
F_pretension_per_screw = F_pretension/N_screws
sigma_per_screw_pretension = screw_lox.comp_tension(F_pretension_per_screw, 'M4')
print(f'{sigma_per_screw_pretension/1e6=}')

# loading per screw given laod with safety factor
F_per_screw = (F_fluid*safety_fluid + F_pretension)/N_screws
sigma_per_screw = screw_lox.comp_tension(F_per_screw, 'M4')
print(f'{sigma_per_screw/1e6=}')

# screw strength given loading per screw and screw factor
F_per_screw_w_safety = F_per_screw*safety_screw

sigma_per_screw_w_safety = screw_lox.comp_tension(F_per_screw_w_safety, 'M4')
print(f'{sigma_per_screw_w_safety/1e6=}')

screw_class = screw_lox.find_strength_class(sigma_per_screw_w_safety)
print(f'{screw_class=}')

safety_screw_actual = screw_lox.yield_strength[screw_class]*1e6/screw_lox.comp_tension(F_per_screw, 'M4')
print(safety_screw_actual)

# safety_screw_actual = screw_lox.yield_strength['9.8']*1e6/screw_lox.comp_tension(F_per_screw, 'M4')
# print(safety_screw_actual)

mu_g = 0.18
mu_k = 0.18
torque = screw_lox.mounting_torque(F_pretension_per_screw, 'M4', mu_g, mu_k)
print(f'{torque=}')