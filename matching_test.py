"""
Test routines for matching functions

created on: 25/07/17

@author: eduardo
"""
import unittest
import numpy as np
import numpy.testing as nptest
import matching


class TestMatching(unittest.TestCase):

    def test_weightNorm(self):
        a = np.array([  ])
        m = 4
        res = matching.weightNorm(a, m, [1,1])
        print(res)

    

if __name__ == '__main__':
    unittest.main()
