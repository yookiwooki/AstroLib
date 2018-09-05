#!/usr/bin/env python3

"""
rkf45.py - Runge-Kutta-Fehlberg 4(5) algorithm
"""

import numpy as np
import numpy.linalg as LA
import pudb

class OptionsRKF(object):
    """ Options for RKF integrator

    tol (float): error tolerance
    h0 (float):  initial step size
    """
    def __init__(self,tol,h0):
        self.tol = tol
        self.h0 = h0

class ResultRKF(object):
    """ Output data from RKF integrator
    t: (numpy.array): m x 1 floats with times corresponding to rows of y
    x: (numpy.array): m x n floats with integrated result

    """
    def __init__(self,t,x):
        self.t = t
        self.x = x

def rkf45(func, tspan, x0, options):
    """ RKF 4(5)th order variable step size numerical integration

    Ref:
    Fehlberg, E. "Low order classical runge-kutta formulas with stepsize
    control and their application to some heat transfer problems." (1969).

    Args:

    Returns:
    ResultRKF: Output data from RKF integrator
    """

    # Set coefficients
    alpha = np.array([2/9, 1/3, 3/4, 1, 5/6])

    beta = np.zeros([5,5])
    beta[0,0] =   np.array([2/9])
    beta[1,0:2] = np.array([1/12, 1/4])
    beta[2,0:3] = np.array([69/128, -243/128, 135/64])
    beta[3,0:4] = np.array([-17/12, 27/4, -27/5, 16/15])
    beta[4,0:5] = np.array([65/432, -5/16, 13/16, 4/27, 5/144])

    c = np.array([1/9, 0, 9/20, 16/45, 1/12])

    c_hat = np.array([47/450, 0, 12/25, 32/225, 1/30, 6/25])

    # Initialization
    f = np.zeros([6,x0.size])
    t = tspan[0]
    x = x0
    h = options.h0

    result = ResultRKF(np.array([tspan[0]]), x0)

    while (t < tspan[1]):

        pudb.set_trace()

        # Evaluate derivatives
        f[0,:] = func(t, x)
        f[1,:] = func(t + alpha[0]*h, beta[0,0]*f[0,:])
        f[2,:] = func(t + alpha[1]*h, beta[1,0]*f[0,:] + beta[1,1]*f[1,:])
        f[3,:] = func(t + alpha[2]*h, beta[2,0]*f[0,:] + beta[2,1]*f[1,:] \
                + beta[2,2]*f[2,:])
        f[4,:] = func(t + alpha[3]*h, beta[3,0]*f[0,:] + beta[3,1]*f[1,:] \
                + beta[3,2]*f[2,:] + beta[3,3]*f[3,:])
        f[5,:] = func(t + alpha[4]*h, beta[4,0]*f[0,:] + beta[4,1]*f[1,:] \
                + beta[4,2]*f[2,:] + beta[4,3]*f[3,:] + beta[4,4]*f[4,:])

        # Evaluate integrated states
        x_hat = x
        for k in range(0,6):
            x_hat = x_hat + h*c_hat[k]*f[k]
        for k in range(0,5):
            x = x + h*c[k]*f[k]

        # Check accuracy
        error = LA.norm(abs(x-x_hat))

        # Save current step
        result.t = np.vstack([result.t, t])
        result.x = np.vstack([result.x, x])

        # Step forward
        t = t + h

    result.x = result.x.transpose()

    return result
