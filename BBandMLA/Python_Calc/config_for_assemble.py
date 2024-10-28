# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:55:30 2024

@author: Ludwig P. Lind
"""

def configuration(E_full, E, x_s, y_s, size, j, setup):
    """
    Parameters
    ----------
    E_full : np.array of shape (2,n_x_full,n_y_full)
        full array of E-field
    E : np.array of shape (2, size, size)
        E-field of single beam for given parameters
    x_s : integer
        start index for adding E onto E_full in central position
    y_s : integer
        same as x_s but in y direction
    size : integer
        size of single E-filed
    j : integer
        index that matches shift p
    setup : string
        string that gives type of setup. options are:
            'cross', 'single', 'X', '3x3', '1x3', '3x1', '5x5'

    Returns
    -------
    E_full: with added E fields of the different beams

    """
    x_e = x_s + size
    y_e = y_s + size
    
    if setup == 'single':
        E_full[:,x_s:x_s+size, y_s:y_s+size] += E
    
    elif setup == 'cross':
        E_full[:,x_s:x_e, y_s:y_s+size] += E
        
        E_full[:,x_s+j:x_e+j, y_s:y_e] += E
        E_full[:,x_s-j:x_e-j, y_s:y_e] += E
        E_full[:,x_s:x_e, y_s+j:y_e+j] += E
        E_full[:,x_s:x_e, y_s-j:y_e-j] += E
    
    elif setup == 'X':
        E_full[:,x_s:x_s+size, y_s:y_s+size] += E
        
        E_full[:,x_s+j:x_e+j, y_s+j:y_e+j] += E
        E_full[:,x_s+j:x_e+j, y_s-j:y_e-j] += E
        E_full[:,x_s-j:x_e-j, y_s+j:y_e+j] += E
        E_full[:,x_s-j:x_e-j, y_s-j:y_e-j] += E
        
    elif setup == '3x3':
        E_full[:,x_s:x_e, y_s:y_s+size] += E
        
        E_full[:,x_s+j:x_e+j, y_s:y_e] += E
        E_full[:,x_s-j:x_e-j, y_s:y_e] += E
        E_full[:,x_s:x_e, y_s+j:y_e+j] += E
        E_full[:,x_s:x_e, y_s-j:y_e-j] += E
        
        E_full[:,x_s+j:x_e+j, y_s+j:y_e+j] += E
        E_full[:,x_s+j:x_e+j, y_s-j:y_e-j] += E
        E_full[:,x_s-j:x_e-j, y_s+j:y_e+j] += E
        E_full[:,x_s-j:x_e-j, y_s-j:y_e-j] += E
        
    elif setup == '1x3':
        E_full[:,x_s:x_e, y_s:y_s+size] += E
        
        E_full[:,x_s:x_e, y_s+j:y_e+j] += E
        E_full[:,x_s:x_e, y_s-j:y_e-j] += E
        
    elif setup == '3x1':
        E_full[:,x_s:x_e, y_s:y_s+size] += E
        
        E_full[:,x_s+j:x_e+j, y_s:y_e] += E
        E_full[:,x_s-j:x_e-j, y_s:y_e] += E
        
    elif setup == '5x5':
        E_full[:,x_s:x_e, y_s:y_s+size] += E
        
        E_full[:,x_s+j:x_e+j, y_s:y_e] += E
        E_full[:,x_s-j:x_e-j, y_s:y_e] += E
        E_full[:,x_s:x_e, y_s+j:y_e+j] += E
        E_full[:,x_s:x_e, y_s-j:y_e-j] += E
        
        E_full[:,x_s+j:x_e+j, y_s+j:y_e+j] += E
        E_full[:,x_s+j:x_e+j, y_s-j:y_e-j] += E
        E_full[:,x_s-j:x_e-j, y_s+j:y_e+j] += E
        E_full[:,x_s-j:x_e-j, y_s-j:y_e-j] += E   
        
        E_full[:,x_s+2*j:x_e+2*j, y_s:y_e] += E
        E_full[:,x_s-2*j:x_e-2*j, y_s:y_e] += E
        E_full[:,x_s:x_e, y_s+2*j:y_e+2*j] += E
        E_full[:,x_s:x_e, y_s-2*j:y_e-2*j] += E
        
        E_full[:,x_s+2*j:x_e+2*j, y_s+2*j:y_e+2*j] += E
        E_full[:,x_s+2*j:x_e+2*j, y_s-2*j:y_e-2*j] += E
        E_full[:,x_s-2*j:x_e-2*j, y_s+2*j:y_e+2*j] += E
        E_full[:,x_s-2*j:x_e-2*j, y_s-2*j:y_e-2*j] += E  
        
        E_full[:,x_s+1*j:x_e+1*j, y_s+2*j:y_e+2*j] += E
        E_full[:,x_s+1*j:x_e+1*j, y_s-2*j:y_e-2*j] += E
        E_full[:,x_s-1*j:x_e-1*j, y_s+2*j:y_e+2*j] += E
        E_full[:,x_s-1*j:x_e-1*j, y_s-2*j:y_e-2*j] += E 
        
        E_full[:,x_s+2*j:x_e+2*j, y_s+1*j:y_e+1*j] += E
        E_full[:,x_s+2*j:x_e+2*j, y_s-1*j:y_e-1*j] += E
        E_full[:,x_s-2*j:x_e-2*j, y_s+1*j:y_e+1*j] += E
        E_full[:,x_s-2*j:x_e-2*j, y_s-1*j:y_e-1*j] += E  

    else:
        print('Invalid Setup Input')
    
    return(E_full)
    