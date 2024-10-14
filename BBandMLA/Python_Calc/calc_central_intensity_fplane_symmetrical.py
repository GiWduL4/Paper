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
from mpmath import nsum, inf, fac, hyp1f1, mpf

# Defining pathes for importing own modules
import sys
sys.path.insert(1, "C:/Users/Ludwig/OneDrive/Dokumente/GitKraken/ATOMICS-python-modules/basic-modules/")
# sys.path.insert(1, "//brain43/groups/Atomics/ATOMICS-python-modules/basic-modules")

import graphical_analysis as ga



"""
B0 and B2; case Z = 0
"""
# Define the series term with additional arguments a, b, z
def B0_term(k, x3, r2):
    return hyp1f1(k+1, 1/2, x3) /fac(k) * (-r2)**k

# Now, use mpmath's nsum to sum the series with additional parameters
def B0_series(rho0, r2):
    # Use nsum to sum from k = 0 to infinity
    return 2 * nsum(lambda k: B0_term(k, -rho0**2, r2), [0, inf])


# Define the series term with additional arguments a, b, z
def B2_term(k, x3, r2):
    return hyp1f1(k+2, 3/2, x3) /fac(k) * (-r2)**k

# Now, use mpmath's nsum to sum the series with additional parameters
def B2_series(rho0, r2):
    # Use nsum to sum from k = 0 to infinity
    return 4 * rho0 * nsum(lambda k: B2_term(k, -rho0**2, r2), [0, inf])


def E_field(params):
    x, y, rho0 = params
    r2 = mpf(x)**2 + mpf(y)**2
    B0 = B0_series(rho0, r2)
    B2 = B2_series(rho0, r2)
    E = np.array([B0 + B2 * x + 1j * B2 * y, B2 * y + 1j * B0 - 1j * B2 * x])*1/np.sqrt(2)
    return E

def E_field_sym(params):
    r, rho0 = params
    r2 = r**2
    B0 = B0_series(rho0, r2)*1/np.sqrt(2)
    # B2 = B2_series(rho0, r2)*1/np.sqrt(2)
    E = [B0, 1j*B0]
    return np.array(E)

def E_field_sym_simplified(params):
    r, rho0 = params
    r2 = r**2
    B0 = B0_series(rho0, r2)
    # B2 = B2_series(rho0, r2)*1/np.sqrt(2)
    E = B0
    return np.array(E)

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

E0 = 1       

n = 300
l = 300

Z = 0

p_list = np.linspace(0,10,l)
rho0_list = np.linspace(0.5,5,n)

# I0 = intensity(E_field((0,0,rho0)))
I = np.zeros((len(rho0_list),len(p_list)))
step = 0



for i,rho0 in enumerate(rho0_list):    
    for j, p in enumerate(p_list):
        E = E_field_sym_simplified((0,rho0)) + 4* E_field_sym_simplified((p,rho0)) #+ E_field((-p,0,rho0)) + E_field((0,p,rho0)) + E_field((0,-p,rho0)) 
#     + E_field(p,-p) + E_field(-p,-p) + E_field(-p,p) + E_field(p,p)
#     + E_field(2*p,0) + E_field(-2*p,0) + E_field(0,2*p) + E_field(0,-2*p)
#     + E_field(2*p,p) + E_field(-2*p,-p) + E_field(2*p,-p) + E_field(-2*p,p)
#     + E_field(p,2*p) + E_field(-p,-2*p) + E_field(-p,2*p) + E_field(p,-2*p))
        I[i,j] = np.abs(E)**2#intensity(Ecross)
    progress = (i+1)/len(rho0_list)
    if progress >= step:
        print('Progress: ' + str(round(progress*100,3)) + ' %')
        step += 0.1
    

"""
Plot
"""

# ga.reel_2D(p_list, rho0_list, I0, xlabel='pitch', ylabel=r'$\rho_0$')
ga.reel_2D(p_list, rho0_list, I, xlabel='pitch', ylabel=r'$\rho_0$', vmax = 10)

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