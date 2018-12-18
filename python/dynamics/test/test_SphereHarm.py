#!/usr/bin/env python3

"""
test_SphereHarm.py - unit tests for SphereHarm.py
"""


import pudb
import unittest
import numpy as np
from dynamics.accel_spharm import SphereHarm

class TestSphereHarm(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
        self.field = 'GRGM1200A'

    def test_SphereHarm(self):
        sph_harm = SphereHarm(self.field)

        #print(sph_harm.sh)
