"""
Test routines for matching functions

created on: 25/07/17

@author: eduardo
"""
import unittest
import numpy as np
import numpy.testing as nptest
import errors

# Initial parameters
ps = 0.003
pm = 0.003
pg = 0.003
eta = 0.01
a0 = 2.0
a1 = 1/80.
theta = .24
protocol = "GHZ"
surface = "toric"

# Create object
errors = errors.Generator(surface=surface, ps=ps, pm=pm, pg=pg, eta=eta,
                          a0=a0, a1=a1, theta=theta, protocol=protocol)

e = errors.get_errors(10, "star", False)
print(e)
