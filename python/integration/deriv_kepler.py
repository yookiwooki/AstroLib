#!/usr/bin/env python3

"""
deriv_kepler.py - Equations of motion for Keplerian dynamics
"""

import pudb
import numpy as np
from astromath.overloads import astro_2norm


def deriv_kepler(t, x):
    """ EOM for two body dynamics

    Args:
    t (float) : time
    x (numpy.array) : 6x1 floats with position and velocity
    Returns:
    f (numpy.array) : 6x1 floats with velocity and acceleration
    """
    # Set mu to 1 for now (TODO: make this configurable)
    mu = 1

    # Initalize output
    f = np.zeros(6)

    # Carry over velocity
    f[0:3] = x[3:6]

    # Calculate and store acceleration
    r_mag = astro_2norm(x[0:3])
    f[3:6] = -(mu/(r_mag**3))*x[0:3]

    return f
