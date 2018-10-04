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
    alpha = np.array([1/4, 3/8, 12/13, 1, 1/2])

    beta = np.zeros([5,5])
    beta[0,0] =   np.array([1/4])
    beta[1,0:2] = np.array([3/32, 9/32])
    beta[2,0:3] = np.array([1932/2197, -7200/2197, 7296/2197])
    beta[3,0:4] = np.array([439/216, -8, 3680/513, -845/4104])
    beta[4,0:5] = np.array([-8/27, 2, -3544/2565, 1859/4104, -11/40])

    c = np.array([25/216, 0, 1408/2565, 2197/4104, -1/5])

    c_hat = np.array([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])

    # Initialization
    f = np.zeros([6,x0.size])
    t = tspan[0]
    x = x0
    h = options.h0

    result = ResultRKF(np.array([tspan[0]]), x0)

    while (t < tspan[1]):

        # Evaluate derivatives
        f[0,:] = h*func(t, x)
        f[1,:] = h*func(t + alpha[0]*h, x + beta[0,0]*f[0,:])
        f[2,:] = h*func(t + alpha[1]*h, x + beta[1,0]*f[0,:] + beta[1,1]*f[1,:])
        f[3,:] = h*func(t + alpha[2]*h, x + beta[2,0]*f[0,:] + beta[2,1]*f[1,:] \
                + beta[2,2]*f[2,:])
        f[4,:] = h*func(t + alpha[3]*h, x + beta[3,0]*f[0,:] + beta[3,1]*f[1,:] \
                + beta[3,2]*f[2,:] + beta[3,3]*f[3,:])
        f[5,:] = h*func(t + alpha[4]*h, x + beta[4,0]*f[0,:] + beta[4,1]*f[1,:] \
                + beta[4,2]*f[2,:] + beta[4,3]*f[3,:] + beta[4,4]*f[4,:])

        # Evaluate integrated states
        x_hat_add = 0
        x_add = 0

        for k in range(0,6):
            x_hat_add = x_hat_add + c_hat[k]*f[k,:]
        for k in range(0,5):
            x_add = x_add + c[k]*f[k,:]

        x_hat = x + x_hat_add
        x = x + x_add

        # Check accuracy
        error = LA.norm(abs(x-x_hat))

        # Save current step
        result.t = np.vstack([result.t, t])
        result.x = np.vstack([result.x, x])

        # Step forward
        t = t + h

    result.x = result.x.transpose()

    return result
