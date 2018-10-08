#!/usr/bin/env python3

"""
test_AstroTime.py - unit tests for AstroTime.py
"""


import pudb
import unittest
import numpy as np
from astrotime.AstroTime import AstroTime

class TestAstroTime(unittest.TestCase):
    def setUp(self):
        self.longMessage = True

    def test_AstroTime(self):
        test_time = AstroTime(2456537.0000000)
        gps_time = test_time.utc2gps(2456537.0000000)
        self.assertEqual(gps_time[0], 1756)
        self.assertEqual(gps_time[1], 43235)

        utc_time = test_time.gps2utc(gps_time)
        self.assertEqual(utc_time, 2456537.0)
