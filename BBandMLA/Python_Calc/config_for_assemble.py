# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:55:30 2024

@author: Ludwig P. Lind
"""

def configuration(E_full, E, x_s, y_s, size, j, setup):
    x_e, y_e = x_s + size, y_s + size

    # Define helper functions for each setup type
    def single_beam():
        E_full[:, x_s:x_e, y_s:y_e] += E

    def cross_beams():
        positions = [(0, 0), (j, 0), (-j, 0), (0, j), (0, -j)]
        for dx, dy in positions:
            E_full[:, x_s+dx:x_e+dx, y_s+dy:y_e+dy] += E

    def x_beams():
        positions = [(0, 0), (j, j), (j, -j), (-j, j), (-j, -j)]
        for dx, dy in positions:
            E_full[:, x_s+dx:x_e+dx, y_s+dy:y_e+dy] += E

    def grid_3x3():
        shifts = [0, j, -j]
        for dx in shifts:
            for dy in shifts:
                E_full[:, x_s+dx:x_e+dx, y_s+dy:y_e+dy] += E

    def grid_5x5():
        shifts = [0, j, -j, 2*j, -2*j]
        for dx in shifts:
            for dy in shifts:
                E_full[:, x_s+dx:x_e+dx, y_s+dy:y_e+dy] += E

    def grid_8x8():
        shifts = [0, j, -j, 2*j, -2*j, 3*j, -3*j, 4*j, -4*j]
        for dx in shifts:
            for dy in shifts:
                if abs(dx) <= 4*j and abs(dy) <= 4*j:  # Ensuring the range for 8x8 grid
                    E_full[:, x_s+dx:x_e+dx, y_s+dy:y_e+dy] += E

    # Map setup names to their corresponding functions
    setups = {
        'single': single_beam,
        'cross': cross_beams,
        'X': x_beams,
        '3x3': grid_3x3,
        '5x5': grid_5x5,
        '8x8': grid_8x8,
    }

    # Execute the chosen setup
    if setup in setups:
        setups[setup]()
    else:
        raise ValueError("Invalid Setup Input")

    return E_full

    