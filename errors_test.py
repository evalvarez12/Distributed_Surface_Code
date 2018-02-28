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
ps = 0.009
pm = 0.009
pg = 0.009
eta = 0.0
a0 = 0.0
a1 = 0.0
theta = 0.0
protocol = "LOCAL_TEST"
surface = "toric"

# Create object
errors = errors.Generator(surface=surface, ps=ps, pm=pm, pg=pg, eta=eta,
                          a0=a0, a1=a1, theta=theta, protocol=protocol)

# symbol = 'IXII_OK'
# print(errors._symbol_to_error(symbol))

e = errors.get_errors(15, "star", False)
print(e[0])
print(e[1])
