# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 18:51:20 2024

@author: Ludwig P. Lind
"""

from mpmath import nsum, inf, fac, hyp1f1
import scipy.special as ss

# Euler's transformation for faster convergence

r = 6.5
r2 = r**2
rho0 = 1.5

x3 = - rho0**2
def alternating_series(k):
    return hyp1f1(k+1, 1/2, x3) /fac(k) * (-r2)**k

sum_result = nsum(alternating_series, [0, inf])
print(sum_result)

A = 0
n = 200
step = 0
for k in range(n):
    A += alternating_series(k)
    progress = (k+1)/n
    if progress >= step:
        print('Progress: ' + str(round(progress*100,3)) + ' %')
        step += 0.1

print(A)