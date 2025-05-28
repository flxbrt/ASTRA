# -*- coding: utf-8 -*-
"""
Created on Tue May 27 21:50:51 2025

@author: felix
"""



from thrustchamber import ThrustChamber


fuel = 'C2H6O'
ox = 'O2'
T_fuel = 400
T_ox = 110
pcc = 20e5
rof = 1
tc = ThrustChamber(fuel, ox, T_fuel, T_ox, pcc)



eta = 1
tc.init_combustion()


tc.set_combustion(p=pcc, rof=rof)
c_star = tc.get_c_star(eta)
print(f'{c_star=:.2f}[m/s] for {rof=} and {pcc/1e5=:.2f}[bar]')

pcc /= 2
tc.set_combustion(p=pcc, rof=rof)
c_star = tc.get_c_star(eta)
print(f'{c_star=:.2f}[m/s] for {rof=} and {pcc/1e5=:.2f}[bar]')