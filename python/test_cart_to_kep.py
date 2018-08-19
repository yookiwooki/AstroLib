#!/usr/bin/env python3

"""
test_cart_to_kep.py - unit tests for kepler.py
"""


import pdb
import unittest
import numpy as np
from numpy import linalg as LA
from cart_to_kep import rv2coe


class TestCartToKep(unittest.TestCase):
    def setUp(self):
        self.mu1 = 1.0 # Canonical (LU^3/sec^2)
        self.mu2 = 398600.4415 # Earth (km^3/sec^2)
        self.tol = 1e-3 # Error tolerance

    def test_noninclined_nonequatorial_elliptical(self):
        rv = np.array([-0.35449, 1.3951 , -0.43398, \
                0.46978, 0.86112, -0.18631])
        coe_ref  = np.array([3.0, 0.815, 2.84563, \
                             2.86945,  .398954, 1.98662])
        coe_test = rv2coe(rv, self.mu1)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for non-inclined, non-equatorial, \
                elliptical orbit (default)')

    def test_elliptical_equatorial(self):
        rv = np.array([1, 0.25, 0, 0.1, 0.5, 0])
        coe_ref  = np.array([0.59513713, 0.78796300, 0.0, \
                             0.0, 0.0, 3.51853886 ])
        coe_test = rv2coe(rv, self.mu1)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for elliptical equatorial orbit')

    def test_circular_inclined(self):
        rv = np.array([-np.sqrt(2), np.sqrt(2), 0, \
            -np.sqrt(2)/4, -np.sqrt(2)/4, 0.5])
        coe_ref  = np.array([2.0, 0.0, 0.785398, \
                             0.0, 2.35619, 0.0])
        coe_test = rv2coe(rv, self.mu1)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for circular inclined orbit')

    def test_circular_equatorial(self):
        rv = np.array([21000*np.sqrt(2), 21000*np.sqrt(2), 0, \
                      -0.05*np.sqrt(self.mu2/210.0), \
                       0.05*np.sqrt(self.mu2/210.0), 0 ])
        coe_ref  = np.array([42000.0, 0, 0, 0, 0, 0.785398])
        coe_test = rv2coe(rv, self.mu2)
        error = LA.norm(coe_ref - coe_test)
        self.assertTrue(error < self.tol,
                'rv2coe fails for circular equatorial orbit')

