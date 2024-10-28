# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:08:47 2024

@author: Ludwig P. Lind
"""

import numpy as np


"""
Setup
"""
E_full = np.zeros((2,401,401), dtype=np.complex128)
x_full = np.linspace(-10,10,401)
y_full = np.linspace(-10,10,401)

"""
Choose parameters
"""
rho0 = 0.924
p = 1
setup = 'cross'

"""
Loading calcs
"""
E_test = np.load('E_single_beam_test.npy')
rho0_list = np.linspace(0.5,5,50)

i = np.argmin(np.abs(rho0_list - rho0))
print('Index ' + str(i) + ' with value ' + r'rho_0 = ' + str(round(rho0_list[i],3)))
 