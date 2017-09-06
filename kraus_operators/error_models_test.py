import unittest
import qutip as qt
import numpy as np
import numpy.testing as nptest
import error_models as errs

# Define variables for tests
psi_test = qt.bell_state("00")
rho_test = psi_test * psi_test.dag()

def Fidelity(rho, psi):
    # Fidelity when one state is a pure state
    return psi.dag() * rho * psi

class TestErrorModels(unittest.TestCase):

    def test_two_qubit_error(self):
        print("----------Test two qubit error-----------")
        p = 1
        rho_noise = errs.two_qubit_gate_noise(rho_test, p, 2, 0, 1)
        print(rho_noise)
        print(Fidelity(rho_noise, psi_test))

    def test_measurement_error(self):
        print("----------Test measurement error-----------")
        p = 1
        psi = qt.snot() * qt.basis(2, 0)
        rho = psi * psi.dag()
        m, rho_coll = errs.measure_single_Xbasis(rho, p)
        print("Measurement: ", m)
        print("Error: p = ", p)
        print(rho_coll)


if __name__ == '__main__':
    unittest.main()
