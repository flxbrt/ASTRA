# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 07:46:25 2025

@author: felix
"""



from CoolProp.CoolProp import PropsSI as psi
    
class continuity:
    def __init__(self, fluid):
        self.fluid = fluid
        
        self.mdot = {'units': ['g/s', 'kg/s'],          # units associated with mass flows
                     'conversion': [1e-3, 1]}
        self.qdot = {'units': ['l/min', 'l/s', 'std l/min', 'm3/h', 'm3/s', 'std cm3/s', 'std cm3/min', 'sccm'], # units associated with voume flows
                     'conversion': [1/1000/60, 1e-3, 1/1000/60, 1/3600, 1, 1e-6, 1e-6/60, 1e-6/60]}
        
        self.mdot = {
           'g/s': 1e-3,
           'kg/s': 1}
        
        self.qdot = {
            'l/min': 1e-3/60,
            'l/s': 1e-3,
            'm3/h': 1/3600,
            'm3/s': 1,
            'std l/min': 1e-3/60,
            'std cm3/s': 1e-6,
            'std cm3/min': 1e-6/60}
        
        self.p_std = 101325             # https://en.wikipedia.org/wiki/International_Standard_Atmosphere
        self.T_std = 288.15             # https://en.wikipedia.org/wiki/International_Standard_Atmosphere
    
    def density(self, p, T):
        self.p = p
        self.T = T
        self.rho = psi('D','P',p,'T',T,self.fluid)
        self.rho_std = psi('D','P',self.p_std,'T',self.T_std,self.fluid)
    
    def to_standard(self, value, unit):
        # converts from any unit into SI units
        # for massflow: kg/s
        # for volume flow: m^3/s
        # does not convert from massflow to volume flow and vice versa
        # self.unit_from = unit
        if unit in self.mdot.keys():
            conversion = self.mdot[unit]
        elif unit in self.qdot.keys():
            conversion = self.qdot[unit]
        else:
            raise ValueError(f'Unit {unit} not defined')
        
        self.standard_value = value*conversion
        
    def to_unit(self, unit, unit_from=None, value=None):
        # converts from SI unit into any unit
        # does not convert from massflow to volume flow and vice versa
        self.unit_to = unit
        if unit_from is not None and value is not None:
            self.to_standard(value, unit_from)
        elif unit_from is not None or value is not None:
            raise ValueError('You have to define unit_from and value')
            
        rho_fac = 1
        if unit in self.mdot.keys():
            conversion = 1/self.mdot[unit]
        elif unit in self.qdot.keys():
            conversion = 1/self.qdot[unit]
            if 'std' in unit:
                rho_fac = self.rho/self.rho_std
        else:
            raise ValueError(f'Unit {unit} not defined')
        self.scaled_value = self.standard_value*conversion*rho_fac
        
    def qdot_mdot(self, p, T, from_to):
        # converts from standard volume into standard mass flow and vice versa
        self.density(p, T)
        
        if from_to == 'q_to_m':
            rho_fac = self.rho
        elif from_to == 'm_to_q':
            rho_fac = 1/self.rho
        else:
            raise ValueError(f'{from_to=} is invalid')
        
        self.standard_value = self.standard_value*rho_fac
        
    def print_all_conversions(self, value, unit, p, T):
        """
        Prints all conversions (massflow + volumeflow) including
        m_dot ↔ q_dot using density(p,T).
        Output is formatted as a clean ASCII table.
        """
    
        # Step 1: save density
        self.density(p, T)
    
        # Step 2: convert input to standard (kg/s or m3/s)
        self.to_standard(value, unit)
        std = self.standard_value
    
        # Determine if input is mass or volume
        input_is_mass = unit in self.mdot
    
        # Step 3: compute both SI standard types
        if input_is_mass:
            mdot_std = std                    # kg/s
            qdot_std = std / self.rho         # m3/s (actual conditions)
        else:
            qdot_std = std                    # m3/s
            mdot_std = std * self.rho         # kg/s
    
        # Step 4: Prepare table rows
        rows = []
        rows.append(["Input", f"{value} {unit}"])
        rows.append(["Medium", f"{self.fluid}"])
        rows.append(["Pressure", f"{self.p/1e5} bar"])
        rows.append(["Temperature", f"{self.T} K"])
        rows.append(["ρ", f"{self.rho:.4g} kg/m³"])
        rows.append(["ρ standard", f"{self.rho_std:.4g} kg/m³"])
        rows.append(["", ""])
    
        rows.append(["--- Mass flow (kg/s basis) ---", ""])
        for u, conv in self.mdot.items():
            val = mdot_std / conv
            rows.append([u, f"{val:.4g}"])
    
        rows.append(["", ""])
        rows.append(["--- Volume flow (m³/s basis) ---", ""])
        for u, conv in self.qdot.items():
            # base conversion
            val = qdot_std / conv
    
            # apply density correction for std units
            if "std" in u:
                val *= (self.rho_std / self.rho)
    
            rows.append([u, f"{val:.4g}"])
    
        # Step 5: Pretty print
        col1_width = max(len(r[0]) for r in rows) + 2
        col2_width = max(len(r[1]) for r in rows) + 2
    
        print("\n" + "="*(col1_width+col2_width))
        print(" ALL UNIT CONVERSIONS ".center(col1_width+col2_width))
        print("="*(col1_width+col2_width))
    
        for r in rows:
            print(f"{r[0]:<{col1_width}} {r[1]:<{col2_width}}")
    
        print("="*(col1_width+col2_width) + "\n")


#%%

hyd = continuity('H2')
p = 40e5
T = 293
value = 3
unit = 'g/s'
hyd.print_all_conversions(value, unit, p, T)

air = continuity('Air')
rof = 200
value *= rof
unit = 'g/s'
air.print_all_conversions(value, unit, p, T)


#%% Option 1 to call this script

value = 3
unit = 'g/s'
fluid = 'N2'
p = 2e5
T = 400

opt1 = continuity(fluid)
opt1.to_standard(value, unit)       # converts the input into SI units --> either into kg/s for a mass flow or into m^3/s for a volume flow
print(f'Input scaled to SI units is {opt1.standard_value}')
opt1.qdot_mdot(p=p, T=T, from_to='m_to_q') # converts the SI unit mass flow in an SI unit volume flow and vice versa
print(f'Input scaled to SI units is {opt1.standard_value}')
opt1.to_unit('std l/min')           # converts from SI unit into anoter unit of the same type (mass flow stays mass flow, volume flow stays volume flow)
print(f'Input scaled to your unit of choice {opt1.scaled_value}')

#%% Option 2 to call the script --> convert in any available unit and print them all

opt2 = continuity(fluid)
opt2.print_all_conversions(value, unit, p, T)

#%% Chat GPT Lösung

# import numpy as np
# from CoolProp.CoolProp import PropsSI as psi

# class Continuity:
#     def __init__(self, fluid):
#         self.fluid = fluid

#         # Mass flows (always to kg/s)
#         self.mdot = {
#             'g/s': 1e-3,
#             'kg/s': 1
#         }

#         # Actual volume flows (always to m^3/s)
#         self.qdot_real = {
#             'l/min': 1e-3/60,
#             'l/s': 1e-3,
#             'm3/h': 1/3600,
#             'm3/s': 1
#         }

#         # Standardized volume flows (to m^3/s at std. conditions)
#         self.qdot_std = {
#             'std l/min': 1e-3/60,
#             'std cm3/s': 1e-6,
#             'std cm3/min': 1e-6/60,
#             'sccm': 1e-6/60
#         }

#         # Standard conditions
#         self.p_std = 1e5
#         self.T_std = 293

#     # -------------------------
#     # Thermodynamic properties
#     # -------------------------
#     def fluid_properties(self, p, T):
#         rho = psi('D', 'P', p, 'T', T, self.fluid)
#         eta = psi('viscosity', 'P', p, 'T', T, self.fluid)
#         sos = psi('speed_of_sound', 'P', p, 'T', T, self.fluid)
#         return rho, eta, sos

#     def compute_density(self, p, T):
#         self.rho = psi('D', 'P', p, 'T', T, self.fluid)
#         self.rho_std = psi('D', 'P', self.p_std, 'T', self.T_std, self.fluid)

#     # -------------------------
#     # Unit conversion
#     # -------------------------
#     def to_standard(self, value, unit):
#         """Convert any unit to internal SI (kg/s or m³/s)."""
#         if unit in self.mdot:
#             factor = self.mdot[unit]
#             self.standard_value = value * factor
#             self.flow_type = "mass"

#         elif unit in self.qdot_real:
#             factor = self.qdot_real[unit]
#             self.standard_value = value * factor
#             self.flow_type = "vol_real"

#         elif unit in self.qdot_std:
#             factor = self.qdot_std[unit]
#             # Convert std-volume → real-volume using density ratio
#             self.standard_value = value * factor
#             self.flow_type = "vol_std"

#         else:
#             raise ValueError(f"Unit {unit} not defined")

#     def to_unit(self, unit, p=None, T=None):
#         """Convert SI (kg/s or m³/s) into any other unit."""
#         rho_fac = 1

#         # Case: converting volume flow (standardized vs real)
#         if unit in self.qdot_std:
#             if p is None or T is None:
#                 raise ValueError("Need p,T for std conversions")
#             self.compute_density(p, T)
#             rho_fac = self.rho / self.rho_std
#             factor = 1 / self.qdot_std[unit]

#         elif unit in self.qdot_real:
#             factor = 1 / self.qdot_real[unit]

#         elif unit in self.mdot:
#             factor = 1 / self.mdot[unit]

#         else:
#             raise ValueError(f"Unit {unit} not defined")

#         return self.standard_value * factor * rho_fac

#     # -------------------------
#     # Universal table output
#     # -------------------------
#     def print_all_conversions(self, value, unit, p, T):
#         self.to_standard(value, unit)
#         self.compute_density(p, T)

#         print("\n------------------------------------------------------------")
#         print(f"   Converted Mass & Volume Flows (Input: {value} {unit})")
#         print("------------------------------------------------------------")
#         print(f"{'Type':12} | {'Unit':12} | {'Value':12}")
#         print("------------------------------------------------------------")

#         # Mass flow outputs
#         if self.flow_type in ["mass"]:
#             for u in self.mdot:
#                 v = self.to_unit(u)
#                 print(f"{'mass':12} | {u:12} | {v:12.4g}")

#         # Volume flow outputs (real)
#         for u in self.qdot_real:
#             v = self.to_unit(u)
#             print(f"{'vol_real':12} | {u:12} | {v:12.4g}")

#         # Volume flow outputs (std)
#         for u in self.qdot_std:
#             v = self.to_unit(u, p, T)
#             print(f"{'vol_std':12} | {u:12} | {v:12.4g}")

#         print("------------------------------------------------------------\n")