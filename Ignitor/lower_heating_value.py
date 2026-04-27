# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:53:26 2026

@author: felix
"""



# my implementation

comb = ct.Solution('gri30.yaml')              # obtain enthalpy at autoignition conditions

# mass averaged temperature

comb.Y = {'H2':1, 'O2':8}

Y_H2 = comb.Y[comb.species_index('H2')]  # mass fraction of H2


comb.TP = 293, 1e5

h_bc = comb.h


# print(comb_ai.Y)

comb.equilibrate('TP')

h_ac = comb.h

print(h_ac - h_bc)

delta_h = h_ac - h_bc   # J/kg mixture


LHV_H2 = abs(delta_h) / Y_H2


# chat gpt

import cantera as ct

# --- Setup ---
gas = ct.Solution('gri30.yaml')

Ta = 293.15      # K
p = 1e5          # Pa

# Stoichiometric H2/O2 mixture
gas.set_equivalence_ratio(phi=1.0, fuel='H2', oxidizer='O2')
gas.TP = Ta, p

# --- Reactants enthalpy ---
h_react = gas.enthalpy_mass   # J/kg (mixture basis)

# Store fuel mass fraction BEFORE reaction
Y_H2 = gas[gas.species_index('H2')].Y[0]

# --- Equilibrate at constant T, p (TP) ---
gas.equilibrate('TP')

# --- Products enthalpy ---
h_prod = gas.enthalpy_mass    # J/kg (mixture basis)

# --- Heat release per kg mixture ---
delta_h_mix = h_prod - h_react   # negative

# --- Convert to per kg of H2 (fuel basis) ---
LHV_H2 = -delta_h_mix / Y_H2     # J/kg fuel

# --- Output ---
print(f"LHV of H2: {LHV_H2/1e6:.2f} MJ/kg")