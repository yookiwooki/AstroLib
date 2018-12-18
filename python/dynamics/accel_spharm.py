#!/usr/bin/env python3

"""
accel_spharm.py - Spherical harmonic gravity acceleration
"""

import pudb
import numpy as np
import requests

class SphereHarm(object):
      """ Specify and calculate spherical harmonics accelerations """

    def __init__(self, field, output=False):

        self.output = output

        # Associate gravity field with filename
        if (field == 'GRGM1200A'):
            fname = 'GRGM1200A.txt'

        # Search data dir for coefficient file
        found = False
        try:
            fh = open('data/' + fname, 'r')
            self.fname = fname
            found = True

        except FileNotFoundError:
            self.download_coeff(fname)
            try:
                fh = open('data/' + fname, 'r')
                self.fname = fname
                found = True
            except FileNotFoundError:
                found = False
                print('ERROR: coefficient file not found locally and'
                      ' and could not be accessed online (SphereHarm)')

        if found:
            self.sh = self.parse(fh)

        fh.close()

    def download_coeff(self, fname):
        """ Download spherical harmonic coefficients from url below """

        if (fname == 'GRGM1200A.txt'):
            url = 'http://pds-geosciences.wustl.edu/grail/' + \
            'grail-l-lgrs-5-rdr-v1/grail_1001/shadr/' + \
            'gggrx_1200a_sha.tab'

        print('SphereHarm.download_coeff(): downloading coefficient'
              ' data from ' + url)
        response = requests.get(url)

        print('SphereHarm.download_coeff(): writing coefficient data'
              ' to ' + fname)
        with open('data/' + fname, 'w') as f:
            f.write(response.text)

    def parse(self, fh):
        """ Read spherical harmonic data from file """

        if self.output:
            print('SphereHarm.parse(): parsing spherical harmonic data')
        header = fh.readline()
        sh = np.zeros([721800,4])

        for idx, line in enumerate(fh):
            line = line.strip()
            columns = line.split(',')
            sh[idx,0] = columns[0]
            sh[idx,1] = columns[1]
            sh[idx,2] = columns[2]
            sh[idx,3] = columns[3]

        return sh

    def sh_accel(self, x, order, degree):
        """ Calculate acceleration due to non-spherical gravity """

