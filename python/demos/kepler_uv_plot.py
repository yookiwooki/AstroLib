#!/usr/bin/env python3

import sys
from os.path import expanduser
import pudb
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

sys.path.insert(0, expanduser("~/AstroLib/python"))

from orbits.kepler_uv import kep_uv, OutputKepUV

"""
kepler_uv_plot.py - demo to show trajectory calculated with kepler_uv
"""

mu = 1
rv0 = np.array([1, 0, 0, 0, 1, 0])
tof = np.linspace(0, 10)

x = np.zeros(tof.size)
y = np.zeros(tof.size)
z = np.zeros(tof.size)

pudb.set_trace()

for i in range(0,tof.size):
    output = kep_uv(rv0, tof[i], mu)
    x[i] = output.rv[0]
    y[i] = output.rv[1]
    z[i] = output.rv[2]

ax = plt.axes(projection='3d')
ax.plot3D(x, y, z)
plt.show()
