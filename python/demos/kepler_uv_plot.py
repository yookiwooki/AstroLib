#!/usr/bin/env python3

import sys
from os.path import expanduser
import pudb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D

# Look for modules in top level of AstroLib
sys.path.insert(0, expanduser("~/AstroLib/python"))

from orbits.kepler_uv import kep_uv, OutputKepUV

"""
kepler_uv_plot.py - demo to show kepler_uv results
"""

mu = 1
rv0 = np.array([1, 0, 0, 0, 1, 0])
tof = np.linspace(0, 5, num=500)

n = tof.size

x = np.zeros(n)
y = np.zeros(n)
z = np.zeros(n)

for i in range(0,n):
    output = kep_uv(rv0, tof[i], mu)
    x[i] = output.rv[0]
    y[i] = output.rv[1]
    z[i] = output.rv[2]


# Plot trajectory
fig1 = plt.figure()

ax11 = fig1.add_subplot(111, projection='3d', \
        xlabel='X Position (LU)', \
        ylabel='Y Position (LU)', \
        zlabel='Z Position (LU)')
traj = ax11.plot(x, y, z, label='Trajectory')
traj_start = ax11.scatter(x[0], y[0], z[0], \
        marker='o', color='k', label='Start')
traj_end = ax11.scatter(x[n-1], y[n-1], z[n-1], \
        marker='x', color='r', label='End')
l = max(np.max(abs(x)),np.max(abs(y)),np.max(abs(z)))*1.1
ax11.set_xlim(-l, l)
ax11.set_ylim(-l, l)
ax11.set_zlim(-l, l)
ax11.legend(loc=2)

# Plot convergence for full tof calculation
fig2 = plt.figure()

ax21 = fig2.add_subplot(211, ylabel='x')
i_total = output.diagnostic['i']
x_store = output.diagnostic['x_store']
x_line = ax21.plot(np.arange(0,i_total+1), x_store, '.-')
ax21.set_xlim(0,i_total)
ax21.xaxis.set_major_locator(MaxNLocator(integer=True))
ax21.grid(True)

ax22 = fig2.add_subplot(212, xlabel='iteration', ylabel='f(x) = TOF - t')
deltat_store = output.diagnostic['deltat_store']
fx_line = ax22.plot(np.arange(0,i_total), deltat_store, 'r.-')
ax22.set_xlim(0,i_total)
ax22.xaxis.set_major_locator(MaxNLocator(integer=True))
ax22.grid(True)

plt.show()

