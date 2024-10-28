# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:08:42 2024

@author: Ludwig P. Lind
"""

import numpy as np
from mpmath import nsum, inf, fac, hyp1f1, mpf


# Defining pathes for importing own modules
import sys
sys.path.insert(1, "C:/Users/Ludwig/OneDrive/Dokumente/GitKraken/ATOMICS-python-modules/basic-modules/")
# sys.path.insert(1, "//brain43/groups/Atomics/ATOMICS-python-modules/basic-modules")

import graphical_analysis as ga

E0 = 1       
Z = 0

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

def E_field2(params):
    x_vals, y_vals, rho0 = params
    r2_vals = [[mpf(x)**2 + mpf(y)**2 for y in y_vals] for x in x_vals]
    B0 = B0_series(rho0, r2_vals)
    B2 = B2_series(rho0, r2_vals)
    E = np.array([
        [[B0[i][j] + B2[i][j] * mpf(x) + 1j * B2[i][j] * mpf(y) for j, y in enumerate(y_vals)] for i, x in enumerate(x_vals)],
        [[B2[i][j] * mpf(y) + 1j * B0[i][j] - 1j * B2[i][j] * mpf(x) for j, y in enumerate(y_vals)] for i, x in enumerate(x_vals)]
    ]) * (1 / np.sqrt(2))
    return np.array(E)

def intensity(Efield):
    return(np.abs(Efield[0])**2 + np.abs(Efield[1])**2)

  

"""
Plot
"""
x = np.linspace(-8,8,201)
y = np.linspace(-8,8,201)
xm, ym = np.meshgrid(x,y)
n = 5
rho0_list = np.linspace(0.5,5,n)
E_list = []

step = 0

for i,rho0 in enumerate(rho0_list):
    E = E_field2((x,y,rho0)) 
    # r2 = xm**2 + ym**2
    # A = np.where(r2>=4.666**2) #computable with numpy
    # E[0][A] = 0
    # E[1][A] = 0
    E_list.append(E)
    Icross = intensity(E)
    ga.reel_2D(x, y, Icross, xlabel='x', ylabel=r'y', vmax = 1)
    progress = (i+1)/len(rho0_list)
    if progress >= step:
        print('Progress: ' + str(round(progress*100,3)) + ' %')
        step += 0.1

np.save('E_single_beam_241028_1633.npy', np.array(E_list))

