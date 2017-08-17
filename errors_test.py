"""
Test routines for matching functions

created on: 25/07/17

@author: eduardo
"""
import unittest
import numpy as np
import numpy.testing as nptest
import errors


cumulativeProbs, errs = errors.processErrors(errors.errTestVec, errors.errorStar1)

class TestMatching(unittest.TestCase):

    def test_weightNorm(self):
        print(cumulativeProbs)
        print(errs[0])
        print(errs[1])
        print(errs[2])

if __name__ == '__main__':
    unittest.main()
