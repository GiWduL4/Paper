# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:08:42 2024

@author: Ludwig P. Lind
"""

import numpy as np
import scipy.integrate as integrate
import math as m
import matplotlib.pyplot as plt
import scipy.special as ss

E0 = 1

def integrand_real(t, Z, rho0, k):
    return(np.real(t*np.exp(-t**2/4*(1+1j*Z))*np.cos(t*rho0)*(t/2)**(2*k)))        
def integrand_imag(t, Z, rho0, k):
    return(np.imag(t*np.exp(-t**2/4*(1+1j*Z))*np.cos(t*rho0)*(t/2)**(2*k)))        

n = 400

Z = 0
k = 30
rho0_list = np.linspace(0,10,n)
"""
Integral
"""
H = np.zeros(n)
K = np.zeros(n)

for j in range(n):
    rho0 = rho0_list[j]
    result0 = integrate.quad(integrand_real, 0, np.inf, args = (Z,rho0,k))
    H[j] = 1/(m.factorial(k)) * result0[0]
    result1 = integrate.quad(integrand_imag, 0, np.inf, args = (Z,rho0,k))
    K[j] = 1/(m.factorial(k)) * result1[0]

"""
Kummer
"""
x3 = -rho0_list**2/(1+1j*Z)
L = 2/(1+1j*Z)**(1+k) * ss.hyp1f1(k+1, 1/2, x3)
M = np.real(L)
N = np.imag(L)

"""
Plot
"""
plt.plot(rho0_list, H)

fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlabel(r'$\rho_0$')
ax.set_ylabel('value')
ax.plot(rho0_list, H, 'c-', label = r'integral real')
ax.plot(rho0_list, K, 'b-', label = r'integral imag')
ax.plot(rho0_list, M, 'm-.', label = r'Kummer real')
ax.plot(rho0_list, N, 'r-.', label = r'Kummer imag')
plt.legend(loc = 'best', prop = {'size':15})
plt.show()