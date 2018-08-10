import unittest
import numpy as np
from faulttolerance.surface_code import SurfaceCode

class Test_SurfaceCode(unittest.TestCase):

    def test_init(self):
        distance = 6
        surface = "toric"
        sc = SurfaceCode(distance=distance,
                         surface=surface)
        print(sc.stars)
        print(sc.plaqs)
        print(sc.tags)
        print(sc.qubits)



if __name__ == "__main__":
    unittest.main()


#1# initialize the surface code object
#1distance = 6
#1surface = "toric"
#1sc = surface_code.SurfaceCode(distance, surface)
#1
#1sc.select_measurement_protocol(0, 0, "single_rounds")
#1
#1
#1print("-------------------Testing random noise------------------------")
#1sc.apply_qubit_error(0.5, 0.5)
#1print(sc.qubits)
#1
#1print("-------------------Testing stabilizer measurements-------------")
#1sc.reset()
#1sc.apply_qubit_error(.3, .3)
#1sc.measure_stabilizer_type("star")
#1sc.measure_stabilizer_type("plaq")
#1print(sc.tags)
#1print(sc.qubits)
#1
#1print("-------------------Testing stabilizer lie-------------")
#1sc.reset()
#1sc._stabilizer_lie("S", .5)
#1sc._stabilizer_lie("P", .5)
#1print(sc.qubits[0])













#
#
#
# class TestSurfaceCode(unittest.TestCase):
#     """Test class for the surface code."""
#
#     def test_structure(self):
#         print("-------------------Testing structure------------------------")
#         stars_position = np.array([[0, 0, 0, 2, 2, 2, 4, 4, 4],
#                                   [1, 3, 5, 1, 3, 5, 1, 3, 5]])
#
#         plaqs_position = np.array([[1, 1, 1, 3, 3, 3, 5, 5, 5],
#                                   [0, 2, 4, 0, 2, 4, 0, 2, 4]])
#
#         np.testing.assert_array_equal(stars_position, sc.stars)
#         np.testing.assert_array_equal(plaqs_position, sc.plaqs)
#
#
#     def test_apply_random_error(self):
#         print("-------------------Testing random noise------------------------")
#         sc._apply_noise_qubit(0.5, 0.5)
#         print(sc.qubits)
#
#     def test_measure_stabilizer(self):
#         print("-------------------Testing stabilizer measurements-------------")
#         sc.reset()
#         sc._apply_noise_qubit(.3, .3)
#         sc.measure_stabilizer_type("star")
#         sc.measure_stabilizer_type("plaq")
#         print(sc.tags)
#         print(sc.qubits)
#
#     def test_stabilizer_lie_all(self):
#         print("-------------------Testing stabilizer lie-------------")
#         sc.reset()
#         sc._stabilizer_lie("S", .5)
#         sc._stabilizer_lie("P", .5)
#         print(sc.qubits[0])
#
#
#     # def test_two_rand_stab_qubits(self):
#     #     print("____________________________")
#     #     test_pos = np.array([[0, 1], [0, 2]])
#     #     a,b = sc._two_rand_stab_qubits(test_pos)
#     #     print(test_pos)
#     #     print(a)
#     #     print(b)
#
#     # def test_
#
#
#
#
# if __name__ == '__main__':
#     unittest.main()
