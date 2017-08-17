"""
Test routines for the surface code implementation

created on: 20/07/17

@author: eduardo
"""
import unittest
import numpy as np
import numpy.testing as nptest
import surface_code

# initialize the surface code object
# do NOT change size, it will ruin the test cases.
size = 3
sc = surface_code.SurfaceCode(size)


class TestSurfaceCode(unittest.TestCase):
    """Test class for the surface code."""

    def test_structure(self):
        starsPosition = np.array([[0, 0, 0, 2, 2, 2, 4, 4, 4],
                                  [0, 2, 4, 0, 2, 4, 0, 2, 4]])

        plaqsPosition = np.array([[1, 1, 1, 3, 3, 3, 5, 5, 5],
                                  [1, 3, 5, 1, 3, 5, 1, 3, 5]])

        np.testing.assert_array_equal(starsPosition, sc.stars)
        np.testing.assert_array_equal(plaqsPosition, sc.plaqs)


    def test_applyRandomError(self):
        print("-------------------Testing random noise------------------------")
        sc._applyNoiseQubit(0.5, 0.5)
        print(sc.qubits)

    def test_measureStabilizer(self):
        print("-------------------Testing stabilizer measurements-------------")
        sc.reset()
        sc._applyNoiseQubit(.3, .3)
        sc.measureStabilizerType("star")
        sc.measureStabilizerType("plaq")
        print(sc.tags)
        print(sc.qubits)

    def test_stabilizerLieAll(self):
        print("-------------------Testing stabilizer lie-------------")
        sc.reset()
        sc._stabilizerLie("S", .5)
        sc._stabilizerLie("P", .5)
        print(sc.qubits[0])


    def test_twoRandStabQubits(self):
        print("____________________________")
        testPos = np.array([[0,1],[0,2]])
        a,b = sc._twoRandStabQubits(testPos)
        print(testPos)
        print(a)
        print(b)





if __name__ == '__main__':
    unittest.main()
