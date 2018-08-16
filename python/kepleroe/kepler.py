#!/usr/bin/env python3

"""
kepler.py - conversion to/from Kepler/Cartesian
"""

import warnings
import numpy as np

K_vec = np.array([0, 0, 1])
local_eps = 1e-9

def rv2coe(rv, mu):
    """ Convert from Cartesian to Keplerian.

    Ref:
    Vallado 4th Ed. pg 113

    Args:
    rv (numpy.array): 6x1 numpy array with pos/vel (LU, LU/TU)
    mu (float):       gravitational constant (LU^3/TU^2)

    Returns:
    numpy.array:      6x1 numpy array with Keplerian orbital elements
                          1: semimajor axis (LU)
                          2: eccentricity
                          3: inclination (radians)
                          4: argument of periapsis (radians)
                          5: right ascension of asc. node (radians)
                          6: true anomaly (radians)
    """
    r_vec = rv[0:3]
    v_vec = rv[3:6]
    r_mag = local_2norm(r_vec)
    v_mag = local_2norm(v_vec)

    # Angular momentum
    h_vec = np.cross(r_vec, v_vec)
    h_mag = local_2norm(h_vec)

    # Node vector
    n_vec = np.cross(K_vec, h_vec)  # Node vector
    n_mag = local_2norm(n_vec)

    # Eccentricity vector
    e_vec = ((v_mag**2 - mu/r_mag)*r_vec - np.dot(r_vec, v_vec)*v_vec)/mu
    e_mag = local_2norm(e_vec)

    # Specific energy
    ksi = v_mag**2/2 - mu/r_mag

    # Semimajor axis and semiparameter
    if abs(1 - e_mag) > local_eps:  # Non-paraolic case
        sma = -mu/(2*ksi)
    else:  # Near parabolic case
        warnings.warn("Eccentricity near zero, setting sma to inf")
        sma = np.inf

    # Inclination
    inc = np.arccos(h_vec[2]/h_mag)

    # Right ascension of the ascending node
    # If J component of n vector is negative, pi<raan<2*pi
    raan = np.arccos(n_vec[0]/n_mag)
    if (n_vec[0] < local_eps):
        raan = 2*np.pi - raan

    # Argument of periapsis
    # If K component of e vector is negative, pi<argp<2*pi
    argp = np.arccos(np.dot(n_vec, e_vec)/(n_mag*e_mag))
    if (e_vec[2] < local_eps):
        argp = 2*np.pi - argp

    # True anomaly
    # Quadrant check: if flight path angle is negative, pi<ta<2*pi
    ta = np.arccos(np.dot(e_vec,r_vec)/(e_mag*r_mag))
    if (np.dot(r_vec, v_vec) < local_eps):
        ta = 2*np.pi - ta

    return np.array([sma, e_mag, inc, argp, raan, ta])


def local_2norm(v_in):
    """ Calculate 2 norm manually.

    This function is included to ensure support for complex numbers,
    which might not work properly in numpy's norm function
    """
    return np.sqrt(v_in[0]**2 + v_in[1]**2 + v_in[2]**2)
