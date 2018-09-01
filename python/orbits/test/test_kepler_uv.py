#!/usr/bin/env python3

"""
test_kepler_uv.py - unit tests for universal varialbes kepler solver
"""

import os
import pickle
import unittest
import pdb
import numpy as np
import numpy.linalg as LA
from orbits.kepler_uv import kep_uv, OutputKepUV

# Define location of test_kepler_uv.py's reference files
modpath = os.path.join(os.path.dirname(__file__), 'ref', '')

def ref_compare(output_test, filename):
    tol = 1e-3

    with open(filename, 'rb') as input:
        output_ref = pickle.load(input)

    if (LA.norm(output_test.rv - output_ref.rv) > tol):
        return False

    for key in output_ref.out_detail:
        resid = output_ref.out_detail[key] - output_test.out_detail[key]
        if (abs(resid) > tol):
            return False

    return True

class TestKeplerUV(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
        self.mu = 1
        self.tof = 5
        self.tol = 1e-3

    def test_case1(self):
        rv0 = np.array([1, 0, 0, 0, 1, 0])
        output_test= kep_uv(rv0, self.tof, self.mu)
        test_pass = ref_compare(output_test,
                modpath + 'kepler_uv_a.pkl')
        self.assertTrue(test_pass)

    def test_case2(self):
        rv0 = np.array([1, 0, 0, 1, 1, 0])
        output_test= kep_uv(rv0, self.tof, self.mu)
        test_pass = ref_compare(output_test,
                modpath + 'kepler_uv_b.pkl')
        self.assertTrue(test_pass)

    def test_case3(self):
        rv0 = np.array([1, 0, 0, 1.414, 0, 0])
        output_test= kep_uv(rv0, self.tof, self.mu)
        test_pass = ref_compare(output_test,
                modpath + 'kepler_uv_c.pkl')
        self.assertTrue(test_pass)

    def test_case4(self):
        rv0 = np.array([1, 0, 0, 0, 2, 0])
        output_test= kep_uv(rv0, self.tof, self.mu)
        test_pass = ref_compare(output_test,
                modpath + 'kepler_uv_d.pkl')
        self.assertTrue(test_pass)


