#!/usr/bin/env python3

"""
AstroTime.py - general time representation and conversion for AstroLib
"""

import pudb
import os.path
import requests
import numpy as np
import time

class AstroTime(object):

    def __init__(self, utc_time):
        self.leapsec = self.get_leapsec(utc_time)

    def utc2gps(self, utc_time):
        """ Convert from UTC time to GPS time

        Ref: Misra and Enge, GPS Signals Measurements and Performance.
        2nd Ed. pg 114

        Args:
        utc_time (float): UTC time in Julian Date format

        Returns:
        np.array: 2x1 floats with GPS week and GPS second
        """
        gps_time = np.zeros([2,1])

        # UTC time in sec since Midnight Jan 5, 1980
        t = (utc_time - 2444244.5)*86400.0 + self.leapsec

        gps_time[0] = int(t/604800)
        gps_time[1] = t % 604800
        return gps_time

    def gps2utc(self, gps_time):
        """ Convert from GPS time to UTC time

        Ref: Misra and Enge, GPS Signals Measurements and Performance.
        2nd Ed. pg 114

        Args:
        gps_time (np.array): 2x1 floats with GPS week and GPS second

        Returns:
        float: UTC time in Julian Date format
        """
        t = gps_time[0]*604800.0 + gps_time[1]

        utc_time = (t - self.leapsec)/86400.0 + 2444244.5
        return utc_time

    def get_leapsec(self, utc_time):
        """ Return the integer number of leap seconds at given time

        Args:
        utc_time (float): UTC time in Julian Date format

        Returns:
        float: Number of leap seconds

        Note:
        A table of leap seconds is stored in the data directory.  This
        function will attempt to open the data file.  If the data file
        is not present, the function will attempt to load the data file
        from the URL below.  If the data file is present, the function
        will check if it is more than six months old, and if so, it will
        download the most up to date data file.
        """

        fname = 'data/leapsec.txt'
        url = 'https://www.ietf.org/timezones/data/leap-seconds.list'
        found = False
        current = False

        print('AstTime.get_leapsec(): Looking for leap second data in ' +
        fname)
        try:
            fh = open(fname,'r')
            found = True

            fmod_time = os.path.getmtime('data/leapsec.txt')
            curr_time = time.time()
            if (curr_time - fmod_time < 6*30*7*24*60*60):
                current = True
            else:
                fh.close()
                print('AstTime.get_leapsec(): leap second data more than six' +
                ' months old')

        except FileNotFoundError:
            print('AstTime.get_leapsec(): leap second data not found')

        if not found or not current:
            print('AstTime.get_leapsec(): Downloading leap second ' +
            'data from ' + url)
            response = requests.get(url)

            print('AstTime.get_leapsec(): Writing leap second data file')
            with open(fname,'w') as f:
                f.write(response.text)
            fh = open(fname,'r')

        print('AstTime.get_leapsec(): Parsing data from ' + fname)

        store_time = 0
        store_leapsec = 10
        for line in fh:
            if not line[0] == '#':
                columns = line.split()

                # Convert time stamp on this line to JD
                line_time = float(columns[0])/86400 + 2415020.5

                if (utc_time > store_time) and (utc_time < line_time):
                    fh.close()
                    return store_leapsec
                store_time = line_time
                store_leapsec = float(columns[1])

        fh.close()
        return float(columns[1])
