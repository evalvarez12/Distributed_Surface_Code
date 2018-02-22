"""
Test routines for the functions in error_models.py

author: Eduardo Villasenor
created-on: 05/06/17
"""
import unittest
import qutip as qt
import error_models as errs

# Define variables for tests
psi_test = qt.bell_state("00")
rho_test = psi_test * psi_test.dag()


class TestErrorModels(unittest.TestCase):

    def test_two_qubit_error(self):
        print("----------Test two qubit error-----------")
        p = .09
        rho_noise = errs.two_qubit_gate_noise(rho_test, p, 2, 0, 1)
        print(rho_noise)
        print(rho_noise.tr())
        print(qt.fidelity(rho_noise, psi_test))

    def test_measurement_error(self):
        print("----------Test measurement error-----------")
        p = .5
        psi = qt.snot() * qt.basis(2, 0)
        rho = psi * psi.dag()
        m, rho_coll = errs.measure_single_Xbasis_random(rho, p)
        print("Measurement: ", m)
        print("Error: p = ", p)
        print(rho_coll)

    def test_env_error(self):
        print("----------Test env dephasing-----------")
        n_steps = 10
        rho_noise1 = errs.env_error(rho_test, 20, .3, 25e-6, 2, [0, 1])
        rho_noise2 = errs.env_error(rho_test, 0, .3, 25e-6, 2, [0, 1])
        # print(rho_noise1)
        # print(rho_noise.tr())
        print(qt.fidelity(rho_noise1, psi_test))
        print(qt.fidelity(rho_noise2, psi_test))



if __name__ == '__main__':
    unittest.main()
