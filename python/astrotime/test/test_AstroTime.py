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

    def test_AstroTime1(self):
        # September 1, 2013 12:00:00 UTC
        ref_utc_time = 2456537.0
        ref_gps_time = np.array([1756, 43235.0])

        self.compare_times(ref_utc_time, ref_gps_time)

    def test_AstroTime2(self):
        # October 8, 2018 18:58:15 UTC
        ref_utc_time = 2458400.2904514
        ref_gps_time = np.array([2022, 154732.0 ])

        self.compare_times(ref_utc_time, ref_gps_time)

    def test_AstroTime3(self):
        # November 5, 1980 03:30:00 UTC
        ref_utc_time = 2444548.6458333
        ref_gps_time = np.array([43, 271819])

        self.compare_times(ref_utc_time, ref_gps_time)

    def compare_times(self, ref_utc_time, ref_gps_time):

        # Instantiating class with utc time
        t1 = AstroTime(ref_utc_time, 'utc')
        gps_time1 = t1.utc2gps(ref_utc_time)

        self.assertAlmostEqual(float(gps_time1[0]), float(ref_gps_time[0]),2)
        self.assertAlmostEqual(float(gps_time1[1]), float(ref_gps_time[1]),2)

        utc_time1 = t1.gps2utc(gps_time1)

        self.assertAlmostEqual(utc_time1, ref_utc_time, 2)

        # Instantiating class with gps time
        t2 = AstroTime(ref_gps_time, 'gps')
        utc_time2 = t2.gps2utc(ref_gps_time)

        self.assertAlmostEqual(utc_time2, ref_utc_time, 2)

        gps_time2 = t2.utc2gps(utc_time2)

        self.assertAlmostEqual(float(gps_time2[0]), float(ref_gps_time[0]),2)
        self.assertAlmostEqual(float(gps_time2[1]), float(ref_gps_time[1]),2)

