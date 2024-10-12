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

def intensity(E_real, E_imag):
    return(np.abs(E_real)**2 + np.abs(E_imag)**2)
    
def E_field_setup(params):
    x, y, rho0 = params
    E = E_field((x,y,rho0)) #+ E_field((x+p,y,rho0)) + E_field((x-p,y,rho0))
    return E

# Vectorized version of electric field calculation for arrays of (x, y)
def E_field_vectorized(X, Y, rho0):
    r2 = X**2 + Y**2  # This is now vectorized, applied to arrays

    B0_vals = np.zeros(X.shape, dtype=complex)
    B2_vals = np.zeros(X.shape, dtype=complex)
    step = 0
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            # Calculate B0 and B2 using mpmath for each (x, y)
            r2_ij = mpf(r2[i, j])
            B0_vals[i, j] = B0_series(rho0, r2_ij)
            B2_vals[i, j] = B2_series(rho0, r2_ij)
        progress = (i+1)/X.shape[0]
        if progress >= step:
            print('Progress: ' + str(round(progress*100,3)) + ' %')
            step += 0.1
    # Now calculate E[0] and E[1] based on B0 and B2
    E_real = B0_vals + B2_vals * X + 1j * B2_vals * Y
    E_imag = B2_vals * Y + 1j * B0_vals - 1j * B2_vals * X

    return E_real, E_imag

# Helper function for parallelization
def compute_intensity(params):
    E = E_field(params)
    return intensity(E)

"""
Optimized Grid Calculation
"""
# Parameters
rho0 = 1.5
p = 1.3  # Using mpmath's arbitrary precision float
nx, ny = 50, 50

# Generate grid
X = np.linspace(-2, 2, nx)
Y = np.linspace(-2, 2, ny)

# Create meshgrid for X and Y
Xm, Ym = np.meshgrid(X, Y)

# Calculate electric field using the vectorized approach
E_r, E_i = E_field_vectorized(Xm, Ym, rho0)

I = intensity(E_r, E_i)

# Create intensity map
# I = np.zeros((nx, ny))
# step = 0
# for i, x in enumerate(X):
#     for j, y in enumerate(Y):
#         E = E_field((x, y, rho0))
#         I[i, j] = intensity(E[0], E[1])
#     progress = (i+1)/nx
#     if progress >= step:
#         print('Progress: ' + str(round(progress*100,3)) + ' %')
#         step += 0.1


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