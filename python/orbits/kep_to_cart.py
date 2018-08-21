#!/usr/bin/env python3

"""
kep_to_cart.py - conversion from Keplerian to Cartesian parameterization
"""

import numpy as np
from astromath.rotations import rot1, rot2, rot3

_local_eps = 1e-9

def coe2rv(coe, mu):
    """ Convert from Keplerian to Cartesian.

    Ref:
    Vallado 4th Ed. pg 118

    Args:
    coe (numpy.array): 6x1 numpy array with Keplerian orbital elements
                           1: semimajor axis (LU)
                           2: eccentricity
                           3: inclination (radians)
                           4: argument of periapsis (radians)
                           5: right ascension of asc. node (radians)
                           6: fast variable (radians)
    mu (float):        gravitational constant (LU^3/TU^2)

    Returns:
    numpy.array:       6x1 numpy array with Cartesian position/velocity
    """
    # Check if circular equatorial
    if ((np.abs(coe[1]) < _local_eps) and (np.abs(coe[2]) < _local_eps)):
        coe[3] = 0
        coe[4] = 0
    # Check if circular inclined
    elif (np.abs(coe[1]) < _local_eps):
        coe[3] = 0
    # Check if elliptical equatorial
    elif (np.abs(coe[2]) < _local_eps):
        coe[4] = 0

    e = coe[1]
    nu = coe[5]
    p = coe[0]*(1 - e**2) # Semiparameter

    # Calculate r and v in perifocal frame
    r_pqw = np.array([p*np.cos(nu)/(1 + e*np.cos(nu)), \
                      p*np.sin(nu)/(1 + e*np.cos(nu)), 0])
    v_pqw = np.array([-np.sqrt(mu/p)*np.sin(nu), \
                       np.sqrt(mu/p)*(e + np.cos(nu)), 0])

    # Calculate total rotation matrix
    r_a = rot3(-coe[3])
    r_b = rot1(-coe[2])
    r_c = rot3(-coe[4])
    r_total = np.matmul(r_c, np.matmul(r_b, r_a))

    r = np.matmul(r_total, r_pqw)
    v = np.matmul(r_total, v_pqw)

    return np.concatenate([r, v])

