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
from mpmath import nsum, inf, fac, hyp1f1, mp, mpf
from multiprocessing import Pool


# Defining pathes for importing own modules
import sys
sys.path.insert(1, "C:/Users/Ludwig/OneDrive/Dokumente/GitKraken/ATOMICS-python-modules/basic-modules/")
# sys.path.insert(1, "//brain43/groups/Atomics/ATOMICS-python-modules/basic-modules")

import graphical_analysis as ga

E0 = 1       

n = 500
k_max = 120

Z = 0
# p_list = [1]#np.linspace(0,5.5,500)
# rho0 = [1.5]#np.linspace(0.5,5,n)
"""
B0 and B2; case Z = 0
"""
# def B0_series(k):
#     return hyp1f1(k+1, 1/2, x3) /fac(k) * (-r2)**k

# sum_result = nsum(alternating_series, [0, inf])
# print(sum_result)


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
    r2 = x**2 + y**2
    B0 = B0_series(rho0, r2)
    B2 = B2_series(rho0, r2)
    E = [B0 + B2 * x + 1j * B2 * y, B2 * y + 1j * B0 - 1j * B2 * x]
    return np.array(E)

def intensity(Efield):
    return(np.abs(Efield[0])**2 + np.abs(Efield[1])**2)
    
def E_field_setup(params):
    x, y, rho0 = params
    E = E_field((x,y,rho0)) #+ E_field((x+p,y,rho0)) + E_field((x-p,y,rho0))
    return E

# Helper function for parallelization
def compute_intensity(params):
    E = E_field(params)
    return intensity(E)

"""
Optimized Grid Calculation
"""
# Parameters
rho0 = mp(1.5)
p = mp(1.3)  # Using mpmath's arbitrary precision float
nx, ny = 20, 20

# Generate grid
X = np.linspace(-3, 3, nx)
Y = np.linspace(-3, 3, ny)

# Create intensity map
I = np.zeros((nx, ny))
step = 0
for i, x in enumerate(X):
    for j, y in enumerate(Y):
        E = E_field((x, y, rho0))
        I[i, j] = intensity(E)
    progress = (i+1)/nx
    if progress >= step:
        print('Progress: ' + str(round(progress*100,3)) + ' %')
        step += 0.1


# Create a list of all grid points as input for parallel computation
# grid_points = [(mpf(x), mpf(y), rho0) for x in X for y in Y]

# # Use multiprocessing to compute intensities in parallel
# with Pool() as pool:
#     I_flat = pool.map(compute_intensity, grid_points)

# # Reshape the flat list back to a 2D grid
# I = np.reshape(I_flat, (nx, ny))

ga.reel_2D(X, Y, I, xlabel='x', ylabel=r'y', vmax = 1)

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