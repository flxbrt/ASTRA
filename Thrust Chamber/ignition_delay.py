# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 09:22:52 2025

@author: felix
"""



import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


def rof_to_phi(rof, rof_stoich):
    # computes the equivalence ratio from the rof
    
    m_ox_rel = 1+1/rof
    m_ox_stoich_rel = 1+1/rof_stoich
    
    lam = m_ox_rel/m_ox_stoich_rel
    
    phi = 1/lam
    
    return phi

# def compute_ignition_delay(mech, fuel, oxidizer, phi, T0, P0):
#     gas = ct.Solution(mech)
#     # Mischungsverhältnis setzen basierend auf stöchiometrische Zusammensetzung
#     gas.TP = T0, P0
#     gas.set_equivalence_ratio(phi, fuel, oxidizer, basis='mole')
#     r = ct.IdealGasReactor(gas)
#     sim = ct.ReactorNet([r])
#     times, temps = [], []
#     while sim.time < 0.01:  # bis max 10 ms
#         sim.step()
#         times.append(sim.time)
#         temps.append(r.T)
#     temps = np.array(temps)
#     t_ign = times[np.argmax(np.gradient(temps))]
#     return t_ign#, times, temps

def ignition_delay(states, species):
    i_ign = states(species).Y.argmax()
    return states.t[i_ign]

def compute_ignition_delay(mech, fuel, oxidizer, phi, T0, P0):
    gas = ct.Solution(mech)
    # Mischungsverhältnis setzen basierend auf stöchiometrische Zusammensetzung
    gas.TP = T0, P0
    gas.set_equivalence_ratio(phi, fuel, oxidizer, basis='mass')
    r = ct.Reactor(contents=gas)
    sim = ct.ReactorNet([r])
    times, temps = [], []
    time_history = ct.SolutionArray(gas, extra=['t'])
    while sim.time < 0.01:  # bis max 10 ms
        t = sim.step()
        times.append(sim.time)
        temps.append(r.T)
        time_history.append(r.thermo.state, t=t)
    temps = np.array(temps)
    t_ign_1 = times[np.argmax(np.gradient(temps))]
    t_ign_2 = ignition_delay(time_history, 'oh')
    return t_ign_1, t_ign_2

def plot_for_validation(tau, T0, pressure, phi):
    """
    Plots ignition delay (tau) vs. 1000/T for constant pressure and phi.

    Parameters:
    - tau: 3D array (shape: (n_T, n_p, n_phi)), but only n_p=1, n_phi=1 here!
    - T0: 1D array of temperatures [K]
    - pressure: scalar, pressure value (for title/legend)
    - phi: scalar, equivalence ratio (for title/legend)
    - unit_pressure: str, unit of pressure (default: Pa)
    - unit_time: str, unit of time (default: s)
    """
    # Only use first pressure and phi
    tau_plot = tau[:, 0, 0]
    inv_T = 1000 / T0

    plt.figure(figsize=(8,6))
    plt.plot(inv_T, tau_plot*1000, 'o-', label=f'p = {pressure/ct.one_atm:.1f} bar\nφ = {phi:.2f}')
    plt.xlabel('1000/T [K⁻¹]', fontsize=12)
    plt.ylabel('Ignition delay τ [ms]', fontsize=12)
    plt.title(f'Ignition Delay vs. Inverse Temperature\n(p = {pressure/ct.one_atm:.1f} bar, φ = {phi:.2f})', fontsize=14)
    plt.yscale('log')
    plt.grid(True, which='both', axis='both', alpha=0.3)
    plt.legend()
    plt.show()


#%%

mech = 'gri30_marinov.yaml'

fuel = 'C2H5OH:1'

oxidizer = 'O2:0.21, N2:0.78'

T0 = np.linspace(1110, 1670, 8)
p0 = np.linspace(3*ct.one_atm, 3*ct.one_atm, 1)
phi0 = np.array([0.5])

tau_1 = np.zeros((len(T0), len(p0), len(phi0)))

'''
cantera schlägt vor ignition delay nach tau_2 zu berechnen https://cantera.org/dev/examples/python/reactors/NonIdealShockTube.html
tau_1 ist von chat_gpt implementiert worden
'''
tau_2 = np.zeros((len(T0), len(p0), len(phi0)))

counter = 0
for index in np.ndindex(tau_1.shape):
    print(counter)
    T = T0[index[0]]
    p = p0[index[1]]
    phi = phi0[index[2]]
    tau_1[index], tau_2[index] = compute_ignition_delay(mech, fuel, oxidizer, phi, T, p)
    counter+=1
    # print(index)
    
plot_for_validation(
    tau_2,
    T0,
    p0[0],  # Nur den ersten Druckwert nehmen
    phi0[0])


# for T in T0:
#     for p in p0:
#         for phi in phi0:
#             pass
            # tau, times, temps = compute_ignition_delay(mech, fuel, oxidizer, phi, T0, P0)



# plotte t_ign über 1000/T als log




