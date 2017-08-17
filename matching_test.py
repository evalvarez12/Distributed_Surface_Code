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

    # def test_proyector(self):

    #     P = proyector_single_qubit(1, 1, 0)
    #     P_ref = qt.Qobj([[0, 1]])
    #     print(P)
    #     print(P_ref)
    #     # self.assertEqual(P, P_ref)
    #
    # def test_proyection(self):
    #     S = qt.bell_state('01')
    #     P = proyector_single_qubit(0, 2, 0)
    #     print(P*S)



if __name__ == '__main__':
    unittest.main()
