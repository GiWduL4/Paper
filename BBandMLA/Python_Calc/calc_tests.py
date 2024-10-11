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

# Defining pathes for importing own modules
import sys
sys.path.insert(1, "C:/Users/Ludwig/OneDrive/Dokumente/GitKraken/ATOMICS-python-modules/basic-modules/")
sys.path.insert(1, "//brain43/groups/Atomics/ATOMICS-python-modules/basic-modules")

import graphical_analysis as ga

E0 = 1       

n = 200
k_max = 60

Z = 0
p_list = np.linspace(0,4.5,200)
rho0_list = np.linspace(0.8,5,n)
"""
B0 and B2; case Z = 0
"""

# first value for central and second for nearest neighbors

# B2 = [0,0]
I0 = np.zeros((len(rho0_list),len(p_list)))
Icross = np.zeros((len(rho0_list),len(p_list)))

for i in range(len(rho0_list)):
    rho0 = rho0_list[i]
    x3 = -rho0**2               
    B00 = 2*ss.hyp1f1(1, 1/2, x3)
    for j in range(len(p_list)):
        p = p_list[j]
        B01 = 0
        for k in range(k_max):
            B01 += 2*ss.hyp1f1(k+1, 1/2, x3)*(-p**2)**k/m.factorial(k)
        I0[i,j] = B00**2
        Icross[i,j] = (B00+4*B01)**2
    # B2[1] += 4*rho0*ss.hyp1f1(k+2, 3/2, x3)*(-p**2)**k/m.factorial(k)

# I0 = B00**2
# Icross = (B00+4*B01)**2

"""
Plot
"""

# ga.reel_2D(p_list, rho0_list, I0, xlabel='pitch', ylabel=r'$\rho_0$')
ga.reel_2D(p_list, rho0_list, Icross, xlabel='pitch', ylabel=r'$\rho_0$', vmax = 0.5)

Imin = np.min(Icross, axis = 0)
rho_opt = rho0_list[np.argmin(Icross, axis = 0)]

fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlabel(r'pitch to waist')
ax.set_ylabel( r'optimal $\rho_0$')

ax.plot(p_list, rho_opt, 'b+', label = r'cross')
# ax.plot(p_list, Imin, 'r-', label = r'cross')
# # ax.plot(p, I0*25, 'g-', label = r'single')
# # ax.plot(rho0_list, N, 'r-.', label = r'Kummer imag')
# plt.legend(loc = 'best', prop = {'size':15})
plt.show()