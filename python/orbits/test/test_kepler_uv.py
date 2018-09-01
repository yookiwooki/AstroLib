#!/usr/bin/env python3

"""
test_kepler_uv.py - unit tests for universal varialbes kepler solver
"""

import unittest
import pdb
import numpy as np
import numpy.linalg as LA
from orbits.kepler_uv import kep_uv, OutputKepUV


class TestKeplerUV(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
        self.mu = 1
        self.tol = 1e-3

    def test_case1(self):
        rv0 = np.array([1, 0, 0, 0, 1, 0])
        tof = 5

        output = kep_uv(rv0, tof, self.mu)

        rv_test = output.rv
        rv_ref = np.array([0.2836621854632262, \
                0.9589242746631377, \
                0.0, \
                -0.9589242746631381, \
                0.2836621854632262, \
                0.0])

        error = LA.norm(rv_test - rv_ref)
        self.assertTrue(error < self.tol,
                'kep_uv fails because of incorrect final rv')
