# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 18:51:20 2024

@author: Ludwig P. Lind
"""

from mpmath import nsum, inf

# Euler's transformation for faster convergence
def alternating_series(k):
    return (-1)**k / (k + 1)

# sum_result = nsum(alternating_series, [0, inf])
# print(sum_result)
A = 0
for k in range(100):
    A += alternating_series(k)

print(A)