import unittest
from operations import *
import qutip as qt
import numpy.testing as nptest

class TestStringMethods(unittest.TestCase):

    def test_tensor_single_operator(self):
        z = qt.sigmaz()
        Z = tensor_single_operator(z, 2, 0)
        # 1 0 0 0
        # 0 1 0 0
        # 0 0 -1 0
        # 0 0 0 -1
        Z_ref = qt.tensor(qt.sigmaz(), qt.qeye(2))
        # print(Z)
        # print(Z_ref)
        self.assertEqual(Z, Z_ref)

    def test_proyector(self):
        print("----------Test proyector-----------")
        P = proyector_single_qubit_Zbasis(1, 1, 0)
        P_ref = qt.Qobj([[0, 1]])
        print(P)
        print(P_ref)
        # self.assertEqual(P, P_ref)

    def test_proyection(self):
        print("----------Test proyection----------")
        S = qt.bell_state('01')
        P = proyector_single_qubit_Zbasis(0, 2, 0)
        print(P*S)

    def test_measurement(self):
        print("----------Test measurement----------")
        S = qt.bell_state('01')
        print(S)
        S = S * S.dag()
        measurement, collapsed_state = measure_single_Zbasis(S, 2, 0, True)
        print(measurement)
        print(collapsed_state)


    # def test_isupper(self):
    #     self.assertTrue('FOO'.isupper())
    #     self.assertFalse('Foo'.isupper())
    #
    # def test_split(self):
    #     s = 'hello world'
    #     self.assertEqual(s.split(), ['hello', 'world'])
    #     # check that s.split fails when the separator is not a string
    #     with self.assertRaises(TypeError):
    #         s.split(2)

if __name__ == '__main__':
    unittest.main()
