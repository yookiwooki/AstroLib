#!/usr/bin/env python3

import numpy as np
from scipy.special import factorial

"""
stumpff.py - Stumpff functions used for universal variables Kepler solver
"""

local_eps = 1.0e-6

def stumpff_c(z):
    """ Calculate Stumpff C variable """

    # Elliptical
    if ((z > 0) and (abs(z) > local_eps)):
        return (1 - np.cos(np.sqrt(z)))/z

    # Parabolic
    elif ((abs(z) <= local_eps)):
        return 0.5 - z/factorial(4) + z**2/factorial(6) \
                - z**3/factorial(8)

    # Hyperbolic
    elif ((z < 0) and (abs(z) > local_eps)):
        return (1 - np.cosh(np.sqrt(-z)))/z

def stumpff_s(z):
    """ Calculate Stumpff S variable"""

    # Elliptical
    if ((z > 0) and (abs(z) > local_eps)):
        return (np.sqrt(z) - np.sin(np.sqrt(z)))/np.sqrt(z**3)

    # Parabolic
    elif ((abs(z) <= local_eps)):
        return 1.0/factorial(3) - z/factorial(5) \
                + z**2/factorial(7) - z**3/factorial(9)

    # Hyperbolic
    elif ((z < 0) and (abs(z) > local_eps)):
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/np.sqrt(-z**3)


