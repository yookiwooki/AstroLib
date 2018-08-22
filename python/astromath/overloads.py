#!/usr/bin/env python3

"""
overloads.py - customized fundamental math functions
"""


import numpy as np


_local_eps = 1e-9

def astro_arccos(ang):
    """ Calculate arccos

    This function includes a guard for cases where floating point error
    leaves the input slightly outside the valid range for arccos.
    """
    if (abs(1 - ang) < _local_eps):
        return 0
    if (abs(1 + ang) < _local_eps):
        return np.pi
    else:
        return np.arccos(ang)

def astro_2norm(v_in):
    """ Calculate 2 norm manually.

    This function is included to ensure support for complex numbers,
    which might not work properly in numpy's norm function
    """
    return np.sqrt(v_in[0]**2 + v_in[1]**2 + v_in[2]**2)

