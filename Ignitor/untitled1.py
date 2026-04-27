# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 10:06:02 2026

@author: felix
"""



from volume_massflow_conversion import continuity


example = continuity('Air')

example.to_standard(1, 'l/s')

print(example.standard_value)

example.qdot_mdot(p=10e5, T=293, from_to='q_to_m')

print(example.rho, example.rho_std)

print(example.standard_value)

example.to_unit('g/s')

print(example.scaled_value)




example = continuity('H2')

example.to_standard(1.7, 'g/s')

print(example.standard_value)

example.qdot_mdot(p=3e5, T=300, from_to='m_to_q')

print(example.rho, example.rho_std)

print(example.standard_value)

example.to_unit('l/s')

print(example.scaled_value)

example.to_unit('std l/min')



            # 'l/min': 1e-3/60,
            # 'l/s': 1e-3,
            # 'm3/h': 1/3600,
            # 'm3/s': 1,
            # 'std l/min': 1e-3/60,
            # 'std cm3/s': 1e-6,
            # 'std cm3/min': 1e-6/60}