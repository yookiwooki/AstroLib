#!/usr/bin/env python3

"""
sc_state.py - store/convert spacecraft states
"""

import numpy as np
from kepler import rv2coe

class State:
    """Store the state of the spacecraft.

    The main function of this class is to convert between state
    parameterizations.
    """

    def __init__(self, state_in, state_type, mu):

        self.state = state_in
        self.state_type = state_type
        self.mu = mu

    def kepler(self):
        if self.state_type == "cartesian":
            return rv2coe(self.state, self.mu)
        if self.state_type == "keplerian":
            return self.state_in
