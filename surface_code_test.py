"""
Test routines for the surface code implementation

created on: 20/07/17

@author: eduardo
"""
import unittest
import numpy as np
import numpy.testing as nptest
import surface_code

size = 3
sc = surface_code.SurfaceCode(size)
class TestSurfaceCode(unittest.TestCase):

    def test_measureStabilizer(self):
        # print(sc.stars)
        print("____________________________")
        sc._applyNoiseQubit(.5,0)
        print(sc.qubits[0])
        print("____________________________")
        sc.measureAllStabilizer("star")
        print(sc.qubits[0])
        print(sc.tags)

    def test_stabilizerLieAll(self):
        print("____________________________")
        sc._stabilizerLieAll(1)
        print(sc.qubits)


    def test_twoRandStabQubits(self):
        print("____________________________")
        print(sc.stars)
        testPos = np.array([[0,1],[0,2]])
        a,b = sc._twoRandStabQubits(testPos)
        print(testPos)
        print(a)
        print(b)



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
