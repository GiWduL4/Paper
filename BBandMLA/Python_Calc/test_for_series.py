# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 18:51:20 2024

@author: Ludwig P. Lind
"""

from mpmath import nsum, inf, fac

# Euler's transformation for faster convergence
def alternating_series(k):
    return (-2)**k / fac(k)

sum_result = nsum(alternating_series, [0, inf])
print(sum_result)

A = 0
n = 500000
step = 0
for k in range(n):
    A += alternating_series(k)
    progress = (k+1)/n
    if progress >= step:
        print('Progress: ' + str(round(progress*100,3)) + ' %')
        step += 0.1

print(A)