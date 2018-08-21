#!/usr/vin/env python3

"""
test_rotations.py - unit tests for rotations.py
"""


import pdb
import unittest
import numpy as np
from numpy import linalg as LA

from astromath.rotations import rot1, rot2, rot3


class TestRotations(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
        self.tol = 1e-9

    def test_rot1(self):
        v = np.array([0, 0, 1])
        ang = np.pi/2
        ans_test = np.matmul(rot1(ang), v)
        ans_ref = np.array([0, 1, 0])
        error = LA.norm(ans_test - ans_ref)
        self.assertTrue(error < self.tol, "rot1 fails")

    def test_rot2(self):
        v = np.array([0, 0, 1])
        ang = np.pi/2
        ans_test = np.matmul(rot2(ang), v)
        ans_ref = np.array([-1, 0, 0])
        error = LA.norm(ans_test - ans_ref)
        self.assertTrue(error < self.tol, "rot2 fails")

    def test_rot3(self):
        v = np.array([1, 0, 0])
        ang = np.pi/2
        ans_test = np.matmul(rot3(ang), v)
        ans_ref = np.array([0, -1, 0])
        error = LA.norm(ans_test - ans_ref)
        self.assertTrue(error < self.tol, "rot3 fails")

