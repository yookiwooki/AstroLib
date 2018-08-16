#!/usr/bin/env python3

import numpy as np
from sc_state import State

rv = np.array([1, 0, 0, 0, 0.5, .1])
mu = 1

test_state = State(rv, "cartesian", mu)
test_kep = test_state.kepler()

print(test_kep)
