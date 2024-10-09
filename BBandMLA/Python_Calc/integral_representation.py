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

Z = 0
k = 0
rho0_list = np.linspace(-2,2,200)

H = np.zeros(200)
for j in range(200):
    rho0 = rho0_list[j]
    result = integrate.quad(integrand_real, 0, np.inf, args = (Z,rho0,k))
    H[j] = 1/(m.factorial(k)) * result[0]


plt.plot(rho0_list, H)
                        