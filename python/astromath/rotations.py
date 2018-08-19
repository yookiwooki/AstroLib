#!/usr/bin/env python3

"""
rotations.py - standard 3D rotation matrices

Note:
Input angles are to be specified in radians.
"""


import numpy as np


def rot1(ang):
    """ Rotation matrix about x or 1 axis"""
    return np.array([[1,            0,           0],
                     [0,  np.cos(ang), np.sin(ang)],
                     [0, -np.sin(ang), np.cos(ang)]])

def rot2(ang):
    """ Rotation matrix about y or 2 axis"""
    return np.array([[np.cos(ang), 0, -np.sin(ang)],
                     [          0, 1,            0],
                     [np.sin(ang), 0,  np.cos(ang)]])

def rot3(ang):
    """ Rotation matrix about z or 3 axis"""
    return np.array([[ np.cos(ang), np.sin(ang), 0],
                     [-np.sin(ang), np.cos(ang), 0],
                     [           0,           0, 1]])

