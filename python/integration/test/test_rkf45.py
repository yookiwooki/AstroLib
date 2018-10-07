#!/usr/bin/env python3

"""
test_rkf45.py - unit tests for rkf45.py
"""


import pudb
import unittest
import numpy as np
from integration.rkf45 import rkf45, OptionsRKF, ResultRKF
from integration.deriv_kepler import deriv_kepler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class TestRKF45(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
        self.options = OptionsRKF(1e-10, 0.001, 0.5, 1.01, 0.8)

    def test_noninc_noneq_ellip(self):

        rv0 = np.array([1, 0, 0, 0.1, 1, 0.1])
        tspan = [0, 6*np.pi]
        result = rkf45(deriv_kepler, tspan, rv0, self.options)

        x = result.x[0,:]
        y = result.x[1,:]
        z = result.x[2,:]

        #fig1 = plt.figure()
        #ax = fig1.add_subplot(111, projection='3d')
        #traj = ax.plot(x, y, z)
        #plt.show()

        #plt.plot(result.h)
        #plt.show()
