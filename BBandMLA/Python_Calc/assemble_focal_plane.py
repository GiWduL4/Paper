# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:08:47 2024

@author: Ludwig P. Lind
"""

import numpy as np
import matplotlib as mpl

import config_for_assemble as cfa
# Defining pathes for importing own modules
import sys
sys.path.insert(1, "C:/Users/Ludwig/OneDrive/Dokumente/GitKraken/ATOMICS-python-modules/basic-modules/")
# sys.path.insert(1, "//brain43/groups/Atomics/ATOMICS-python-modules/basic-modules")

import graphical_analysis as ga

# Graphical SetUps
mpl.rc("figure",figsize=(12,9))
mpl.rc("xtick", labelsize = 18)
mpl.rc("ytick", labelsize = 18)
mpl.rc("axes", labelsize = 20)

"""
Functions
"""
def intensity(Efield):
    return(np.abs(Efield[0])**2 + np.abs(Efield[1])**2)


"""
Choose parameters
"""
rho0 = 0.924
p = 8
setup = '3x3'

"""
Loading calcs
"""
H = np.load('test.npz')

x = H['arr_0'] # np.linspace(-5,5,201)
y = H['arr_1'] # np.linspace(-5,5,201)
size = len(x)
dx = x[1]-x[0]
rho0_list = H['arr_2']    #np.linspace(0.5,5,50)
E_test = H['arr_3']       #np.load('E_single_beam_test.npy') #WIP: add x,y,rho0_list to npy

i = np.argmin(np.abs(rho0_list - rho0))
print('Index ' + str(i) + ' with value ' + r'rho_0 = ' + str(round(rho0_list[i],3)))

j = round(p/dx)
print('Pixel number ' + str(j) + ' with distance ' + r'p = ' + str(round(j*dx,3)))

print('Setup: ' + setup)

"""
Setup for E_full
"""
l_max = 20
n = 2*round(l_max/dx)+1
E_full = np.zeros((2,n,n), dtype=np.complex128)
x_full = np.linspace(-l_max,l_max,n)
y_full = np.linspace(-l_max,l_max,n)


"""
Calculating
"""
E = E_test[i]
x_s = int((n-size)/2)
y_s = int((n-size)/2)

E_full = cfa.configuration(E_full, E, x_s, y_s, size, j, setup)

I = intensity(E_full)

"""
Plotting
"""

ga.reel_2D(x_full, y_full, I, xlabel='x/w', ylabel=r'y/w', vmax = 1.54)

 