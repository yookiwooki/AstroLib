#!/usr/bin/env python3

"""
refgen_kepler_uv.py - generates reference data for test_kepler_uv.py
"""

import sys
from os.path import expanduser
import pickle
import numpy as np
import numpy.linalg as LA
# Look for modules in top level of AstroLib
sys.path.insert(0, expanduser("~/AstroLib/python"))
from orbits.kepler_uv import kep_uv, OutputKepUV


def refgen(rv0, filename):
    mu = 1
    tof = 5

    print('Saving to: ' + filename)
    data = kep_uv(rv0, tof, mu)

    with open(filename, 'wb') as output:
        pickle.dump(data, output, pickle.HIGHEST_PROTOCOL)

    del data

    with open(filename, 'rb') as input:
        data_stored = pickle.load(input)

    print(data_stored.rv)
    print(data_stored.out_detail)
    print(data_stored.diagnostic)
    print('\n')

# Case A: Circular Orbit
print('Case A: Circular Orbit')
rv0_a = np.array([1, 0, 0, 0, 1, 0])
refgen(rv0_a, 'ref/kepler_uv_a.pkl')

# Case B: Hyperbolic Orbit 1
print('Case B: Hyperbolic Orbit 1')
rv0_b = np.array([1, 0, 0, 1, 1, 0])
refgen(rv0_b, 'ref/kepler_uv_b.pkl')

# Case C: Rectilinear Orbit
print('Case C: Rectilinear Orbit')
rv0_c = np.array([1, 0, 0, 1.414, 0, 0])
refgen(rv0_c, 'ref/kepler_uv_c.pkl')

# Case D: Hyperbolic Orbit 2
print('Case D: Hyperbolic Orbit 2')
rv0_d = np.array([1, 0, 0, 0, 2, 0])
refgen(rv0_d, 'ref/kepler_uv_d.pkl')

