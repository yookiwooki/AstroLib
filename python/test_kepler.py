#!/usr/bin/env python3

"""
test_kepler.py - unit tests for kepler.py
"""


import unittest
import numpy as np
from numpy import linalg as LA
from kepler import rv2coe


class TestKepler(unittest.TestCase):
    def setUp(self):
        self.mu = 1.0 # Canonical (LU^3/sec^2)
        self.tol = 1e-5 # Error tolerance

    def test_noninclined_nonequatorial_elliptical(self):
        rv = np.array([1, 0, 0, 0, 0.5, .1])
        coe_ref  = np.array([0.57471264, 0.74, 0.19738556, \
                             3.14159265,  0.0, 3.14159265])
        coe_test = rv2coe(rv, self.mu)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for non-inclined, non-equatorial, \
                elliptical orbit (default)')

    def test_elliptical_equatorial(self):
        rv = np.array([1, 0, 0, 0, 0.5, 0])
        coe_ref  = np.array([0.57471264, 0.74, 0.19738556, \
                             3.14159265,  0.0, 3.14159265])
        coe_test = rv2coe(rv, self.mu)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for elliptical equatorial orbit')

    def test_circular_inclined(self):
        rv = np.array([1, 0, 0, 0, 0.5, 0.1])
        coe_ref  = np.array([0.57471264, 0.74, 0.19738556, \
                             3.14159265,  0.0, 3.14159265])
        coe_test = rv2coe(rv, self.mu)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for circular inclined orbit')

    def test_circular_equatorial(self):
        rv = np.array([1, 0, 0, 0, 0.5, .1])
        coe_ref  = np.array([0.57471264, 0.74, 0.19738556, \
                             3.14159265,  0.0, 3.14159265])
        coe_test = rv2coe(rv, self.mu)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for circular equatorial orbit')

