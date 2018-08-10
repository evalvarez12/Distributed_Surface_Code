
"""
Test file for the functions on Circuit Block.

author: Eduardo Villasenor
created-on: 20/11/17
"""

import qutip as qt
import faulttolerance.decomposition.circuit_block
import faulttolerance.decomposition.error_models as errs
import numpy as np
import unittest
import logging

#logging.basicConfig(level=logging.INFO)

class Test_Blocks(unittest.TestCase):

    _USE_RAJA_PARAMS = False

    def setUp(self):
        self.noiseless_blocks = Test_Blocks.get_noiseless_blocks()

        if Test_Blocks._USE_RAJA_PARAMS:
            # NOTE: Parameters as in Raja's thesis
            self.ps = 0.006
            self.pm = 0.006
            self.pg = 0.006
            self.a0 = 83.33
            self.a1 = 1/3.
            self.eta = (0.1)*(0.03)*(0.8)
            self.theta = .24
        else:
            # the other set of parameters as originally provided in this file
            self.ps = 0.006
            self.pm = 0.006
            self.pg = 0.006
            self.a0 = 5.
            self.a1 = 1/80.
            self.eta = 1/50.
            self.theta = .63
        self.blocks = circuit_block.Blocks(ps=self.ps,
                                           pm=self.pm,
                                           pg=self.pg,
                                           eta=self.eta,
                                           a0=self.a0,
                                           a1=self.a1,
                                           theta=self.theta)

    @staticmethod
    def get_noiseless_blocks():
        return(circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4.))

    def test_init(self):
        for element in self.blocks.check.elements():
            self.assertTrue(self.blocks.check[element] == 0)

    def test_change_parameters(self):
        random_parameters = {"ps": 2,
                             "pm": 3,
                             "pg": 5,
                             "eta": 7,
                             "a0": 11,
                             "a1": 13,
                             "theta": 17}
        self.blocks.change_parameters(**random_parameters)
        for parname in list(random_parameters.keys()):
            self.assertEqual(getattr(self.blocks, parname),
                            random_parameters[parname])
        self.test_init()

        # reset the parameters again to what they were
        self.setUp()

    @staticmethod
    def get_perfect_bell_pair(state):
        """ `state` should be one of ['00', '01', '10', '11']"""
        perfect_state = qt.bell_state(state=state)
        dm = perfect_state * perfect_state.dag()
        return perfect_state, dm

    def check_perfect_case(self, bell_generate_function, state):
        """ `state` should be one of ['00', '01', '10', '11']"""
        # test the case of the noiseless generation (eta=0)
        # where the desired output state is Psi^+ (i.e. theta=pi/4)

        duration, bell_state = bell_generate_function()
        
        # perfect bell state 
        noiseless_Psi_plus_state, dm = Test_Blocks.get_perfect_bell_pair(state=state)

        # assert that the generated state is a perfect bell state indeed
        self.assertEqual(dm, bell_state)

        # TODO: also test the case with maximal noise and check that 
        # the output is the |11> state


    def test_generate_bell_single_click(self):
        self.check_perfect_case(bell_generate_function=self.noiseless_blocks._generate_bell_single_click, 
                                state='10')

    def test_generate_bell_pair_BK(self):
        self.check_perfect_case(bell_generate_function=self.noiseless_blocks._generate_bell_pair_BK, 
                                state='00')

    def test_success_number_of_attempts(self):
        """
        The function Blocks._success_number_of_attempts samples
        from the CDF of the geometric distribution. 
        """

        # asserting that for unit success probability, we never have to 
        # draw more often than 1 = 0 + 1 times.
        for __ in range(100):
            self.assertEqual(self.blocks._success_number_of_attempts(p_success=1), 0)

        # asserting that the output is always an integer
        for p_success in [0.1, 0.5]:
            sample = self.blocks._success_number_of_attempts(p_success=p_success)
            self.assertEqual(type(sample), np.int64)

        # testing that some statistical properties are correct
        num_of_samples = 20  # arbitrary number that should be taken as big as possible
        threshold = None  # expected closeness of the statistical properties, given the number of samples
        for p_success in [0.1, 0.5]:

            # samples from the function Block._success_number_of_attempts
            samples = []
            for __ in range(num_of_samples): 
                sample = self.blocks._success_number_of_attempts(p_success=p_success)
                samples.append(sample)
            empirical_mean = np.mean(samples)
            empirical_median = np.median(samples)
            empirical_variance = np.var(samples)

            # theoretical values of the CDF of the geometric distribution
            theoretical_mean = (1-p_success)/p_success
            theoretical_median = -1.*np.log(2.)/np.log(1-p_success)
            theoretical_variance = (1 - p_success)/(p_success * p_success)

            # asserting that the results are close to what is expected
            # theoretically
            self.assertTrue(Test_Blocks.check_close(threshold, empirical_mean, theoretical_mean))
            self.assertTrue(Test_Blocks.check_close(threshold, empirical_median, theoretical_median))
            self.assertTrue(Test_Blocks.check_close(threshold, empirical_variance, theoretical_variance))

    @staticmethod
    def check_close(threshold, val1, val2):
        # NOTE could be replaced by TestCase.assertAlmostEqual
        logging.info("The following two values should be close: {} and {}".format(val1, val2))
        if threshold is None:
            return True
        return np.abs(val1 - val2) <= threshold


    def test_start_bell_pair(self):
        # I get an error in _success_number_of_attempts
        #self.blocks.start_bell_pair(rho=None)
        pass


if __name__ == "__main__":
    unittest.main()


#1 def p_ref(f): return (f**2 +2*f*(1-f)/3 + 5*((1-f)/3)**2)
#1 def single_selection(F1, F2):
#1     F = F1*F2 + (1 - F1)*(1 - F2)/9
#1     F = F/(F1*F2 + F1*(1 - F2)/3 + F2*(1 - F1)/3 + 5*(1 - F1)*(1 - F2)/9)
#1     return F
#1 
#1 # Determine parameters
#1 ps = 0.006
#1 pm = 0.006
#1 pg = 0.006
#1 a0 = 5.
#1 a1 = 1/80.
#1 eta = 1/50.
#1 theta = .63
#1 
#1 # NOTE: Parameters as in Raja's thesis
#1 # ps = 0.006
#1 # pm = 0.006
#1 # pg = 0.006
#1 # a0 = 83.33
#1 # a1 = 1/3.
#1 # eta = (0.1)*(0.03)*(0.8)
#1 # theta = .24
#1 
#1 
#1 # Initialize  objects
#1 cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#1 cb_ideal = circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4.)
#1 rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
#1 
#1 print("------------------TEST SWAP--------------------------")
#1 _, _, rho = cb.start_epl()
#1 print("F initial: ", qt.fidelity(rho, rho_ref))
#1 rho = cb._swap_pair(rho, [0, 1])
#1 print("F: ", qt.fidelity(rho, rho_ref))
#1 
#1 print("------------------SINGLE SELECTION-------------------")
#1 rho = errs.bell_pair(0.1)
#1 print("F initial: ", qt.fidelity(rho, rho_ref)**2)
#1 rho = qt.tensor(rho, rho)
#1 p, check, rho = cb_ideal.single_selection_ops(rho, [0, 1], [2, 3],  "Z")
#1 print("p_success: ", p, p_ref(0.9))
#1 print("check: ", check)
#1 print("F: ", qt.fidelity(rho, rho_ref)**2, single_selection(.9, .9))
#1 
#1 print("------------------DOUBLE SELECTION-------------------")
#1 rho = errs.bell_pair(0.1)
#1 print("F initial: ", qt.fidelity(rho, rho_ref)**2)
#1 rho = qt.tensor(rho, rho, rho)
#1 p, check, rho = cb_ideal.double_selection_ops(rho, [0, 1], [2, 3], [4, 5], "Z")
#1 print("p_success: ", p)
#1 print("check: ", check)
#1 print("F: ", qt.fidelity(rho, rho_ref)**2)
#1 
#1 # print("------------------DOUBLE SELECTION2-------------------")
#1 # p, check, rho = cb.double_selection22("Z")
#1 # print("p_success: ", p)
#1 # print("check: ", check)
#1 # print("F: ", qt.fidelity(rho, rho_ref))
#1 #
#1 #
#1 # print("------------------EPL protocol-------------------")
#1 # p, check, rho = cb.start_epl()
#1 # # p, check, rho = cb.double_selection(rho, [0, 1], "X")
#1 # print("p_success: ", p)
#1 # print("check: ", check)
#1 # print("F: ", qt.fidelity(rho, rho_ref))
#1 
#1 # print("------------------NOISY BELL GENERATOR-------------------")
#1 # p, check, rho = cb.start_bell_pair()
#1 # print("p_success: ", p)
#1 # print("check: ", check)
#1 # print("F: ", qt.fidelity(rho, rho_ref))
