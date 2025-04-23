# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 14:27:45 2024

@author: felix
"""



from CoolProp.CoolProp import PropsSI as psi



T_init = 293
p_init = psi('P','T',T_init,'Q',1,'N2O')
p_end = 30e5

kappa = psi('Cpmass','P',p_init,'Q',1,'N2O')/psi('Cvmass','P',p_init,'Q',1,'N2O')


V_ratio = (p_end/p_init)**(1/kappa)




p_after = p_init*0.5**kappa/1e5



#%%%

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

# Define the gas and temperature
gas = 'NitrousOxide'
temperature = 300  # Temperature in Kelvin

# Define pressure range (from 1 atm to 10 atm in Pascal)
pressures = np.linspace(101325, 60 * 101325, 100)  # Pressure range from 1 atm to 10 atm

# List to store heat capacity ratio values
gamma_values = []

# Calculate Cp/Cv for each pressure value
for pressure in pressures:
    cp = CP.PropsSI('Cpmass', 'T', temperature, 'P', pressure, gas)
    cv = CP.PropsSI('Cvmass', 'T', temperature, 'P', pressure, gas)
    gamma = cp / cv
    gamma_values.append(gamma)

# Plot the heat capacity ratio as a function of pressure
plt.plot(pressures / 101325, gamma_values)
plt.xlabel('Pressure (atm)')
plt.ylabel('Heat Capacity Ratio (γ)')
plt.title(f'Heat Capacity Ratio of {gas} at {temperature} K')
plt.grid(True)
plt.show()
