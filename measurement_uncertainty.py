# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 13:28:36 2024

@author: felix
"""



import numpy as np
import matplotlib.pyplot as plt





# m_dot = 0.48 corresponds to 0.8bar pressure differential

m_dot = 0.4

Kv = 1.8
rho = 800
rho0 = 1000
p1 = 0.8e5
p2 = 0
sigma_p1 = sigma_p2 = np.linspace(0.001e5, 0.1e5, 100)
deltap = p1 - p2

sigma_m_dot_abs = rho*Kv*np.sqrt(rho0/rho/1e5)/2/np.sqrt(p1-p2)*np.sqrt(sigma_p1**2+sigma_p2**2)/3600
sigma_m_dot_dif = rho*Kv*np.sqrt(rho0/rho/1e5)/2/np.sqrt(deltap)*np.sqrt(sigma_p1**2)/3600


plt.plot(sigma_p1/1e5, sigma_m_dot_abs/m_dot*100, label='abs')
plt.plot(sigma_p1/1e5, sigma_m_dot_dif/m_dot*100, label='dif')
plt.legend()
plt.grid()
