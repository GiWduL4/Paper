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
from decimal import Decimal

# Defining pathes for importing own modules
import sys
sys.path.insert(1, "C:/Users/Ludwig/OneDrive/Dokumente/GitKraken/ATOMICS-python-modules/basic-modules/")
# sys.path.insert(1, "//brain43/groups/Atomics/ATOMICS-python-modules/basic-modules")

import graphical_analysis as ga

E0 = 1       

n = 500
k_max = 30

Z = 0
p_list = [1]#np.linspace(0,5.5,500)
rho0 = [1.5]#np.linspace(0.5,5,n)
"""
B0 and B2; case Z = 0
"""
def B0_pre(rho0):
    # Ensure rho0 is an array (even if it's a single scalar value)
    rho0 = np.atleast_1d(rho0)
    
    # Initialize the output as a 2D array (for each rho0, you get an array of size k_max)
    B0_prefactor = np.zeros((len(rho0), k_max))
    
    for i, r in enumerate(rho0):
        x3 = -r**2
        print(x3)
        for k in range(k_max):
            B0_prefactor[i, k] = Decimal(ss.hyp1f1(k+1, 1/2, x3)) / Decimal(m.factorial(k))
    
    # If rho0 was a scalar, return a 1D array instead of 2D
    # if B0_prefactor.shape[0] == 1:
    #     return B0_prefactor[0]
    return B0_prefactor

def B0_calc(r2):
    B0 = 0
    for k in range(k_max):
        B0 += B0pre[:,k]*(-r2)**k 
        # print(type(B0))
    B0 = 2 * B0
    return(B0)

def B2_pre(rho0):
    # Ensure rho0 is an array (even if it's a single scalar value)
    rho0 = np.atleast_1d(rho0)
    
    # Initialize the output as a 2D array (for each rho0, you get an array of size k_max)
    B2_prefactor = np.zeros((len(rho0), k_max))
    
    for i, r in enumerate(rho0):
        x3 = -r**2
        for k in range(k_max):
            B2_prefactor[i, k] = Decimal(ss.hyp1f1(k+2, 3/2, x3)) / Decimal(m.factorial(k))
    
    # If rho0 was a scalar, return a 1D array instead of 2D
    # if B2_prefactor.shape[0] == 1:
    #     return B2_prefactor[0]
    return B2_prefactor

def B2_calc(r2):
    B2 = 0
    for k in range(k_max):
        B2 += B2pre[:,k]*(-r2)**k 
        # print(type(B0))
    B2 = 2 * B2
    return(B2)



def E_field(x,y):
    E = [0,0]
    r2 = x**2 + y**2
    B0 = B0_calc(r2)
    B2 = B2_calc(r2)
    E[0] = B0 + B2*x + 1j*B2*y
    E[1] = B2*y + 1j* B0 -1j*B2*x
    return(np.array(E))

def intensity(Efield):
    return(np.abs(Efield[0])**2 + np.abs(Efield[1])**2)
# first value for central and second for nearest neighbors

# B2 = [0,0]
# I0 = np.zeros((len(rho0_list),len(p_list)))
# Icross = np.zeros((len(rho0_list),len(p_list)))

# for i in range(len(rho0_list)):
#     rho0 = rho0_list[i]
#     x3 = -rho0**2               
#     B00 = 2*ss.hyp1f1(1, 1/2, x3)
#     for j in range(len(p_list)):
#         p = p_list[j]
#         B01 = 0
#         for k in range(k_max):
#             B01 += 2*ss.hyp1f1(k+1, 1/2, x3)*(-p**2)**k/m.factorial(k)
#         I0[i,j] = B00**2
#         Icross[i,j] = (B00+4*B01)**2
    # B2[1] += 4*rho0*ss.hyp1f1(k+2, 3/2, x3)*(-p**2)**k/m.factorial(k)

# I0 = B00**2
# Icross = (B00+4*B01)**2

B0pre = B0_pre(rho0)
B2pre = B2_pre(rho0)
print('prefactors calculated')

# I0 = intensity(E_field(0,0))
# Icross = np.zeros((len(rho0),len(p_list)))
# step = 0
# for j, p in enumerate(p_list):
#     Ecross = E_field(0,0) + E_field(p,0) + E_field(-p,0) + E_field(0,p) + E_field(0,-p)
#     Icross[:,j] = intensity(Ecross)
#     progress = (j+1)/len(p_list)
#     if progress >= step:
#         print('Progress: ' + str(round(progress*100,3)) + ' %')
#         step += 0.1
    

"""
Plot
"""
x = np.linspace(-3,3,400)
y = np.linspace(-3,3,400)
xm, ym = np.meshgrid(x,y)

p = 1.3

I0 = intensity(E_field(xm,ym))
E = E_field(xm,ym)#+ E_field(xm-p,ym) +E_field(xm+p,ym) 
Icross = intensity(E)
ga.reel_2D(x, y, Icross, xlabel='x', ylabel=r'y', vmax = 1)

# ga.reel_2D(p_list, rho0_list, I0, xlabel='pitch', ylabel=r'$\rho_0$')
# ga.reel_2D(p_list, rho0, Icross, xlabel='pitch', ylabel=r'$\rho_0$', vmax = 10)

# Imin = np.min(Icross, axis = 0)
# rho_opt = rho0_list[np.argmin(Icross, axis = 0)]

# fig, ax = plt.subplots()
# plt.grid(True)
# ax.set_xlabel(r'pitch to waist')
# ax.set_ylabel( r'optimal $\rho_0$')

# ax.plot(p_list, rho_opt, 'b+', label = r'cross')
# ax.plot(p_list, Imin, 'r-', label = r'cross')
# # ax.plot(p, I0*25, 'g-', label = r'single')
# # ax.plot(rho0_list, N, 'r-.', label = r'Kummer imag')
# plt.legend(loc = 'best', prop = {'size':15})
# plt.show()