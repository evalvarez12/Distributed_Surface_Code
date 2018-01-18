import qutip as qt
import numpy as np
import error_models as errs
import operations as ops


class EntanglementGeneration():
    """
    Class containing the theoretical functions
    for generating entanglement through optical cavities.
    """

    def __init__(self, p_em, p_ps, p_det, theta):
        """Init function."""
        self.eta = p_em * p_ps * p_det
        self.theta = theta

    def change_parameters(self, p_em, p_ps, p_det, theta):
        """Change class parameters."""
        self.eta = p_em * p_ps * p_det
        self.theta = theta

    def _success_number_of_attempts(self, p_success):
        # Up to 20 tries for success
        i = np.arange(10000000)
        # print(p_success)
        d = self._distribution(p_success, i)
        return np.random.choice(i, 1, p=d)[0]
        # return 0

    def _distribution(self, p, n):
        # Distribution for the probability in the number of tries of
        # each event
        return p*(1-p)**n

    def _bell_pair_click(p, theta):
        s = np.sin(theta)**2
        r = ((1 - p)*s)/(1 - p*s)
        state = qt.bell_state('01') * qt.bell_state('01').dag()
        noise = qt.tensor(qt.basis(2, 1), qt.basis(2, 1))
        noise = noise * noise.dag()
        return (1-r)*state + r*noise

    def _p_success_single_click(self):
        """Probaility of success single click protocol."""
        s = np.sin(self.theta)**2
        c = np.cos(self.theta)**2
        p_success = 2*s*self.eta*(c + s*(1 - self.eta))
        return p_success

    def generate_bell_single_click(self):
        # Probaility of success
        p_success = self._p_success_single_click()

        # This circuit number of steps and time it took
        attempts = self._success_number_of_attempts(p_success) + 1
        time = self.time_lookup["bell_pair"] * attempts

        # Update check
        self.check["bell_pair"] += 1
        # self.check["time"] += time

        # Generate noisy bell pair
        bell = errs._bell_pair_click(self.eta, self.theta)
        return time, bell

    def _p_success_BK(self):
        """Probaility of success Barret-Kok protocol."""
        s = np.sin(self.theta)**2
        r = (1 - self.eta)*s/(1 - self.eta*s)
        p_success = (1 - r)*self.eta**2
        return p_success

    def generate_bell_pair_BK(self):
        # Probaility of success
        p_success = self._p_success_BK()

        # This circuit number of steps
        attempts = self._success_number_of_attempts(p_success) + 1
        time = self.time_lookup["bell_pair"] * attempts

        # Update check
        self.check["bell_pair"] += 1

        # Generate noisy bell pair
        bell = qt.bell_state('00') * qt.bell_state('00').dag()
        return time, bell
