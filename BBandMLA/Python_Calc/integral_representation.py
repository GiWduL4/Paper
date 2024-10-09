# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:08:42 2024

@author: Ludwig P. Lind
"""

import numpy as np
import scipy.integrate as integrate
import math as m
import matplotlib.pyplot as plt

E0 = 1

def integrand_real(t, Z, rho0, k):
    return(np.real(t*np.exp(-t**2/4*(1+1j*Z))*np.cos(t*rho0)*(t/2)**(2*k)))        
def integrand_imag(t, Z, rho0, k):
    return(np.imag(t*np.exp(-t**2/4*(1+1j*Z))*np.cos(t*rho0)*(t/2)**(2*k)))        

Z = 0
k = 1
rho0_list = np.linspace(-2,2,200)

H = np.zeros(200)
K = np.zeros(200)
for j in range(200):
    rho0 = rho0_list[j]
    result0 = integrate.quad(integrand_real, 0, np.inf, args = (Z,rho0,k))
    H[j] = 1/(m.factorial(k)) * result0[0]
    result1 = integrate.quad(integrand_imag, 0, np.inf, args = (Z,rho0,k))
    K[j] = 1/(m.factorial(k)) * result1[0]


plt.plot(rho0_list, H)

fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlabel(r'$\rho_0$')
ax.set_ylabel('value')
ax.plot(rho0_list, H, label = r'integral real')
ax.plot(rho0_list, H, label = r'integral imag')
                        