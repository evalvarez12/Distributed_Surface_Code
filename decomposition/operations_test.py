import unittest
import qutip as qt
import numpy as np
import numpy.testing as nptest
from operations import *

# Define variables for tests
N = 3
pos = 1


def assert_Equivalent(A, B):
    """
    Obtain a global_phase between the two operators.
    """
    # Numpy array representations
    a = A.full()
    b = B.full()

    # Get the first non zero element
    non_zero_ind = np.nonzero(a)
    x = non_zero_ind[0][0]
    y = non_zero_ind[1][0]

    # Find the phase if theres any
    phase = b[x, y] / a[x, y]
    if phase * A == B:
        return True
    return False



class TestOperations(unittest.TestCase):

    def test_tensor_single_operator(self):
        print("----------Test tensor_single_operator-----------")
        z = qt.sigmaz()
        Z = tensor_single_operator(z, N, pos)
        Z_ref = qt.rz(np.pi, N, pos)
        # print(Z)
        # print(Z_ref)
        self.assertTrue(assert_Equivalent(Z, Z_ref))

    def test_proyector_dimRed(self):
        print("----------Test proyector w dimRed-----------")
        P = projector_single_qubit_Zbasis_dimRed(0, N, pos)
        P_ref = qt.Qobj([[1, 0]])
        P_ref = tensor_single_operator(P_ref, N, pos)
        # print(P)
        # print(P_ref)
        self.assertEqual(P, P_ref)

    def test_proyector(self):
        print("----------Test proyector -----------")
        P = projector_single_qubit_Zbasis(1, N, pos)
        P_ref = qt.basis(2, 1) * qt.basis(2, 1).dag()
        P_ref = tensor_single_operator(P_ref, N, pos)
        # print(P)
        # print(P_ref)
        self.assertEqual(P, P_ref)

    def test_measurement_Zbasis(self):
        print("----------Test measurement Z basis----------")
        rho = qt.bell_state('00') * qt.bell_state('00').dag()
        measurement, rho = random_measure_single_Zbasis(rho, 2, 0, True)
        if measurement == 1:
            rho_ref = qt.basis(2, 0) * qt.basis(2, 0).dag()
        if measurement == -1:
            rho_ref = qt.basis(2, 1) * qt.basis(2, 1).dag()
        self.assertEqual(rho, rho_ref)

    def test_measurement_Xbasis(self):
        print("----------Test measurement X basis----------")
        psi = qt.snot() * qt.basis(2, 1)
        rho = qt.tensor(psi, psi)
        rho = rho * rho.dag()
        measurement, rho = random_measure_single_Xbasis(rho, 2, 1, True)
        rho_ref = psi * psi.dag()
        self.assertEqual(rho, rho_ref)


if __name__ == '__main__':
    unittest.main()
