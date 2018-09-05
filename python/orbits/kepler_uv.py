#!/usr/bin/env python3


"""
kepler_uv.py - solve Kepler's problem using universal variables
"""


import pdb
import warnings
import numpy as np
from astromath.overloads import astro_2norm
from orbits.stumpff import stumpff_c, stumpff_s


class OutputKepUV(object):
    """ Output data storage class for kep_uv

    rv (numpy.array):  6x1 floats with final pos/vel (LU, LU/TU)
    out_detail (dict): detailed output data
    diagnostic (dict): detailed diagnostic data
    """
    def __init__(self, rv, out_detail, diagnostic):
        self.rv = rv
        self.out_detail = out_detail
        self.diagnostic = diagnostic

def kep_uv(rv0, tof, mu, rs_tol=1e-12, check_tol=1e-5, imax=100):
    """ Universal variables Kepler solver

    This algorithm predicts the position/velocity of a spacecraft
    with two body dynamics, given the initial position/velocity and
    a time of flight.  The universal variables approach presented by
    Bate, Mueller, and White is particularaly well suited for both
    elliptic, parabolic, and hyperbolic orbits, as well as the edge
    cases between them.

    Ref:
    Fundamentals of Astrodynamics (Bate, Mueller, and White)
    Fundamentals of Ast. and Applications (Vallado)
    Celestial Mechanics II Lecture Notes (Dr. Russell, UT Austin)

    Args:
    rv0 (numpy.array): 6x1 floats with initial pos/vel (LU, LU/TU)
    tof (float):       time of flight (TU)
    mu (float):        gravitational parameter (LU^3/TU^2)
    rs_tol (float):    root solver tolerance
    check_tol (float): post root solver accuracy check on f/g
    imax (int):        root solver max iterations

    Returns
    OutputKepUV:       Output data storage class (see class def.)
    """
    # Step 1: get pos/vel magnitude and SMA
    r0_mag = astro_2norm(rv0[0:3])
    v0_mag = astro_2norm(rv0[3:6])

    a = 0.5*(-mu/(v0_mag**2/2 - mu/r0_mag))
    alpha = -v0_mag**2/mu + 2/r0_mag # 1/a

    # Step 2: Solve UV TOF equation with newton step root solve
    # Initial guess - see Vallado KEPLER algorithm
    if (abs(tof) < 1.0e-6): # Guess zero if TOF near zero
        x = 0
    else:
        if (alpha > 1.0e-6): # Elliptical
            x = np.sqrt(mu)*tof*alpha
        elif (abs(alpha) < 1.0e-6): # Parabolic
            p = (astro_2norm(np.cross(rv0[0:3] ,rv0[3:6])))**2/mu
            s = 1/(np.arctan(3*np.sqrt(mu/(p**3))*tof))/2
            w = np.arctan(np.tan(s)**(1.0/3.0))
            x = np.sqrt(p)*2*1/(np.tan(2*w))
        else: # Hyperbolic
            x = np.sign(tof)*np.sqrt(-a)*np.log(-2*mu*alpha*tof/ \
                    (np.dot(rv0[0:3], rv0[3:6]) + np.sign(tof)* \
                    np.sqrt(-mu*a)*(1 - r0_mag*alpha)))

    deltax = 1
    i = 0

    x_store = []
    z_store = []
    deltat_store = []

    while ((np.abs(deltax) > rs_tol) and (i < imax)):
        z = x**2*alpha
        x_store.append(x)
        z_store.append(z)

        c = stumpff_c(z)
        s = stumpff_s(z)

        r_mag = (x**2)*c + (np.dot(rv0[0:3], rv0[3:6])/np.sqrt(mu))*x* \
                (1 - z*s) + r0_mag*(1 - z*c)

        t = (1/np.sqrt(mu))*((x**3)*s + (np.dot(rv0[0:3], rv0[3:6])/ \
                np.sqrt(mu))*x**2*c + r0_mag*x*(1 - z*s))

        deltat_store.append(tof - t)

        deltax = np.sqrt(mu)*(tof-t)/r_mag
        x = x + deltax

        i = i + 1

    x_store.append(x)
    z_store.append(x**2*alpha)

    if (i == imax):
        warnings.warn("Maximum iterations reached in root solver")

    # Step 3: calculage f, g, f_dot, g_dot
    f = 1 - x**2*c/r0_mag
    g = t - x**3*s/np.sqrt(mu)
    f_dot = np.sqrt(mu)*x*(z*s - 1)/(r0_mag*r_mag)
    g_dot = 1 - x**2*c/r_mag

    # Check accuracy condition here
    if (np.abs(1-(f*g_dot - f_dot*g)) > check_tol):
        warnings.warn("f/g accuracy condition not met.")

    # Step 4: compute final position and velocity
    rv = np.zeros(6)
    rv[0:3] = f*rv0[0:3] + g*rv0[3:6]
    rv[3:6] = f_dot*rv0[0:3] + g_dot*rv0[3:6]

    # Store detailed final data
    out_detail = {}
    out_detail['c'] = c
    out_detail['s'] = s
    out_detail['f'] = f
    out_detail['g'] = g
    out_detail['f_dot'] = f_dot
    out_detail['g_dot'] = g_dot
    out_detail['x'] = x
    out_detail['z'] = z

    # Store detailed diagnostic data
    diagnostic = {}
    diagnostic['i'] = i
    diagnostic['x_store'] = np.asarray(x_store)
    diagnostic['z_store'] = np.asarray(z_store)
    diagnostic['deltat_store'] = np.asarray(deltat_store)

    return OutputKepUV(rv, out_detail, diagnostic)
