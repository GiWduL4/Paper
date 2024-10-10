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

n = 400
k_max = 10

Z = 0
k = 1
p = 2
rho0= np.linspace(0.8,1.1,n)
"""
B0 and B2; case Z = 0
"""
x3 = -rho0**2
# first value for central and second for nearest neighbors
B0 = [0,0]
# B2 = [0,0]
B0[0] = 2*ss.hyp1f1(1, 1/2, x3)
# B2[0] = 3*rho0**ss.hyp1f1(2, 3/2, x3)
for k in range(k_max):
    B0[1] += 2*ss.hyp1f1(k+1, 1/2, x3)*(-p**2)**k/m.factorial(k)
    # B2[1] += 4*rho0*ss.hyp1f1(k+2, 3/2, x3)*(-p**2)**k/m.factorial(k)

I0 = B0[0]**2
Icross = (B0[0]+4*B0[1])**2

"""
Plot
"""
fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlabel(r'$\rho_0$')
ax.set_ylabel('value')

ax.plot(rho0, Icross, 'b-', label = r'cross')
ax.plot(rho0, I0, 'g-', label = r'single')
# ax.plot(rho0_list, N, 'r-.', label = r'Kummer imag')
plt.legend(loc = 'best', prop = {'size':15})
plt.show()