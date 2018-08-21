#!/usr/bin/env python3

"""
test_kep_to_cart.py - unit tests for kep_to_cart.py
"""


import pdb
import unittest
import numpy as np
from numpy import linalg as LA
from orbits.kep_to_cart import coe2rv


class TestKepToCart(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
        self.mu1 = 1.0 # Canonical (LU^3/TU^2)
        self.mu2 = 398600.4415 # Earth (km^3/sec^2)
        self.tol = 1e-2 # Error tolerance

    def test_noninc_noneq_ellip(self):
        rv_ref = np.array([-0.35449, 1.3951 , -0.43398, \
                0.46978, 0.86112, -0.18631])
        coe  = np.array([3.0, 0.815, 2.84563, \
                             2.86945,  .398954, 1.98662])
        rv_test = coe2rv(coe, self.mu1)
        error = LA.norm(rv_ref - rv_test)
        self.assertTrue(error < self.tol,
                'coe2rv fails for non-inclined, non-equatorial, \
                elliptical orbit (default)')

    def test_ellip_eq(self):
        rv_ref = np.array([1, 0.25, 0, 0.1, 0.5, 0])
        coe  = np.array([0.59513713, 0.78796300, 0.0, \
                            3.51853886075583, 0.0,-3.2735602 ])
        rv_test = coe2rv(coe, self.mu1)
        error = LA.norm(rv_ref - rv_test)
        self.assertTrue(error < self.tol,
                'coe2rv fails for elliptical equatorial orbit')

    def test_circ_inc(self):
        rv_ref = np.array([-np.sqrt(2), np.sqrt(2), 0, \
            -np.sqrt(2)/4, -np.sqrt(2)/4, 0.5])
        coe  = np.array([2.0, 0.0, 0.785398, \
                             0.0, 2.35619, 0.0])
        rv_test = coe2rv(coe, self.mu1)
        error = LA.norm(rv_ref - rv_test)
        self.assertTrue(error < self.tol,
                'coe2rv fails for circular inclined orbit')

    def test_circ_eq(self):
        rv_ref = np.array([21000*np.sqrt(2), 21000*np.sqrt(2), 0, \
                      -0.05*np.sqrt(self.mu2/210.0), \
                       0.05*np.sqrt(self.mu2/210.0), 0 ])
        coe  = np.array([42000.0, 0, 0, 0, 0, 0.785398])
        rv_test = coe2rv(coe, self.mu2)
        error = LA.norm(rv_ref - rv_test)
        self.assertTrue(error < self.tol,
                'coe2rv fails for circular equatorial orbit')

