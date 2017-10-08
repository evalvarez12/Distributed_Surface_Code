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
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.0
system_size = 4
function = "LOCAL"
surface = "planar"
# Create object
errors = errors.Generator(surface=surface, ps=ps, pm=pm,
                          pg=pg, pn=pn, protocol=function)

e=errors.get_errors(5, "X", True)
print(e)
