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
            'std m3/s': 1,
            'std m3/h': 1/3600,
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
                # val *= (self.rho_std / self.rho)
                val *= (self.rho / self.rho_std)
    
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
        



if __name__ == '__main__':
    check = continuity('H2')
    check.to_standard(1, 'g/s')
    check.qdot_mdot(1e5, 293, 'm_to_q')
    check.to_unit('std l/min')
    print(check.scaled_value)
    
    check = continuity('Air')
    check.to_standard(30.1, 'g/s')
    check.qdot_mdot(1e5, 293, 'm_to_q')
    check.to_unit('std l/min')
    print(check.scaled_value)
