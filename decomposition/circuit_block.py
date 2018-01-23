"""
Individual circuits involved the remote generation of entanglement
and in different steps of purification protocols.

author: Eduardo Villasenor
created-on: 20/11/17
"""
import qutip as qt
import numpy as np
import collections
import itertools
import operations as ops
import error_models as errs


class Blocks:
    """
    Class for holding all protocols.
    Each circuit block returns the resulting state, number of steps used
    and the probability of success if applicable.

    Parameters
    -----------
    ps - single qubit gate error
    pm - single qubit measurement error
    pg - two-qubit gate error
    pn - network error
    """

    def __init__(self, ps, pm, pg, eta, a0, a1, theta):
        """Init function."""
        # Set the parameters to all faulty opearations
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.eta = eta
        self.a0 = a0
        self.a1 = a1
        self.theta = theta
        # Dictionary to carry the check of all the operations made so far
        self.check = collections.Counter({"bell_pair": 0,
                                          "two_qubit_gate": 0,
                                          "single_qubit_gate": 0,
                                          "measurement": 0,
                                          "time0": 0,
                                          "time1": 0,
                                          "time": 0})

        # Lookup tables for the time it takes to make each operaation
        self.time_lookup = {"bell_pair": 6e-6,
                            "two_qubit_gate": 200e-6,
                            "single_qubit_gate": 200e-6,
                            "measurement": 200e-6}

    def change_parameters(self, ps, pm, pg, eta, a0, a1, theta):
        """Function to change the parameters of all the operations."""
        # Reset all the parameters for the faulty operations
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.eta = eta
        self.a0 = a0
        self.a1 = a1
        self.theta = theta
        self._reset_check()

    def _reset_check(self):
        """Reset to 0 all vaulues from check."""
        # Reset check dictionary
        self.check = collections.Counter(dict.fromkeys(self.check, 0))

    def _generate_bell_single_click(self):
        # Generate a Bell pair using the single click protocol.
        # Only takes into account the error introducced by photon loss.
        # Probaility of success
        s = np.sin(self.theta)**2
        c = np.cos(self.theta)**2
        p_success = 2*s*self.eta*(c + s*(1 - self.eta))

        # This circuit number of steps and time it took
        attempts = self._success_number_of_attempts(p_success) + 1
        time = self.time_lookup["bell_pair"] * attempts

        # Update check
        self.check["bell_pair"] += 1
        # self.check["time"] += time

        # Generate noisy bell pair
        bell = errs.bell_pair_click(self.eta, self.theta)
        return time, bell

    def _generate_bell_pair_BK(self):
        # Generate a Bell pair using the Barret-Kok protocol.
        # Probaility of success
        s = np.sin(self.theta)**2
        r = (1 - self.eta)*s/(1 - self.eta*s)
        p_success = (1 - r)*self.eta**2

        # This circuit number of steps
        attempts = self._success_number_of_attempts(p_success) + 1
        time = self.time_lookup["bell_pair"] * attempts

        # Update check
        self.check["bell_pair"] += 1

        # Generate noisy bell pair
        bell = qt.bell_state('00') * qt.bell_state('00').dag()
        return time, bell

    def _generate_noisy_plus(self):
        # Prepare a noisy |+> state, error is the same
        # as measurement error.
        # This circuit time
        # NOTE: Used time of a single qubit gate
        time = self.time_lookup["single_qubit_gate"]

        # Generate noisy plus
        plus = qt.snot() * qt.basis(2, 0)
        plus = plus * plus.dag()
        minus = qt.snot() * qt.basis(2, 1)
        minus = minus * minus.dag()
        plus = (1 - self.pm) * plus + self.pm * minus

        return time, plus

    def _success_number_of_attempts(self, p_success):
        # Draw a number from the Distribution to simulate
        # the number of attempts required until success.
        # Up to 10000 tries for success
        # NOTE: This is expensive, adjust manually as required
        i = np.arange(1000000)
        # print(p_success)
        d = self._distribution(p_success, i)
        return np.random.choice(i, 1, p=d)[0]
        # return 0

    def _distribution(self, p, n):
        # Distribution for the probability in the number of tries of
        # each event
        return p*(1-p)**n

    def _append_noisy_plus(self, rho):
        # Append noisy |+> state to the total state rho
        # Generate noisy plus
        time, plus = self.generate_noisy_plus()

        # Apply environmental error
        rho = errs.env_error_all(rho, 0, self.a1, time)
        # Noisy plus tate is attached at the end of the complete state
        rho = qt.tensor(rho, plus)
        return rho

    def _append_bell_pair(self, rho):
        # Append single click Bell pair state to the total state rho
        # Generate raw Bell pair to the state.
        time, bell = self._generate_bell_single_click()
        self.check["time"] += time
        self.check["time0"] += time
        # Apply environmental error
        rho = errs.env_error_all(rho, self.a0, self.a1, time)

        # Bell state is attached at the end of the complete state
        rho = qt.tensor(rho, bell)
        return rho

    def _append_epl(self, rho):
        # Append a Bell pair generated using the EPL protocol to the
        # total state rho
        # Copy a backup of the check dictS
        check_backup = self.check.copy()
        # Get EPl and number of attempts
        p_success, _, rho_epl = self.start_epl()
        attempts = self._success_number_of_attempts(p_success) + 1

        # Multiply check elements
        for k in self.check:
            self.check[k] *= (attempts)
        time = self.check["time0"]
        self.check += check_backup

        # Dephase the rest of the state and join them
        rho = errs.env_error_all(rho, self.a0, self.a1, time)
        rho = qt.tensor(rho, rho_epl)
        return rho

    def _get_two_qubit_gates(self, N, controls, targets, sigma):
        # Construct two qubit C-sigma gates.
        gates = []
        if sigma == "X":
            for i in range(len(controls)):
                gates += [qt.cnot(N, controls[i], targets[i])]

        if sigma == "Z":
            for i in range(len(controls)):
                gates += [qt.cphase(np.pi, N, controls[i], targets[i])]

        return gates

    def _collapse_single(self, rho, pos, project, basis):
        # Measure a single qubit in the state.
        # NOTE: Dimension is reduced after collapse
        N = len(rho.dims[0])

        if basis == "X":
            rho = errs.single_qubit_gate_noise(rho, self.ps, N, pos)
            rho = errs.measure_single_Xbasis_forced(rho, self.pm, project, N, pos)
        elif basis == "Z":
            rho = errs.measure_single_Zbasis_forced(rho, self.pm, project, N, pos)
        return rho

    def _swap_noise(self, rho, pos):
        # Apply the noise induced by the one way SWAP gate
        # NOTE: use only two CNOTs to perform a SWAP
        # Swap noise is only single qubit gate because one of the states
        # because a one way Swap gate is used
        N = len(rho.dims[0])

        # Apply noise
        for i in range(2):
            rho = errs.single_qubit_gate_noise(rho, self.pg, N, pos)
        return rho

    def _swap_pair(self, rho, pair):
        # Apply the noise due to SWAP on two states
        # NOTE: use only two CNOTs to perform a SWAP
        self.check["two_qubit_gate"] += 2
        self.check["time"] += self.time_lookup["two_qubit_gate"]*2
        self.check["time1"] += self.time_lookup["two_qubit_gate"]*2
        rho = self._swap_noise(rho, pair[0])
        rho = self._swap_noise(rho, pair[1])
        return rho

    def _collapse_ancillas_X(self, rho, ancillas_pos, projections):
        # Measure the ancillas in the X basis in parallel in each node.
        # NOTE: All ancillas are collapsed in parallel, Hadamard operations
        # are used to measure on X basis
        self.check["measurement"] += 1
        self.check["single_qubit_gate"] += 1

        time = (self.time_lookup["measurement"]
                + self.time_lookup["single_qubit_gate"])
        self.check["time"] += time
        self.check["time1"] += time
        # Apply environmental error
        rho = errs.env_error_all(rho, 0, self.a1, time)

        # Calculate probability of success using the same considerations as
        # in single selection
        N_ancillas = len(ancillas_pos)
        N = len(rho.dims[0])
        if N_ancillas == 2:
            p_success = ops.p_success_single_sel(rho, N, ancillas_pos)
        elif N_ancillas == 4:
            ancillas_pos1 = ancillas_pos[:2]
            ancillas_pos2 = ancillas_pos[2:]
            p_success = ops.p_success_double_sel(rho, N,
                                                 ancillas_pos1,
                                                 ancillas_pos2)
        else:
            p_success = 1

        # Collapse the qubits in parrallel
        # Sort list to be able to reduce dimension and keep track of positions
        ancillas_pos = sorted(ancillas_pos)
        for i in range(N_ancillas):
            pos = ancillas_pos[i] - i
            rho = self._collapse_single(rho, pos,
                                        projections[i], "X")

        return p_success, rho

    def _collapse_ancillas_Z(self, rho, ancillas_pos, projections):
        # Measure the ancillas in the X basis in parallel in each node.
        # NOTE: All ancillas are collapsed in parallel, Hadamard operations
        # are used to measure on X basis
        self.check["measurement"] += 1

        time = self.time_lookup["measurement"]
        self.check["time"] += time
        self.check["time1"] += time
        # Apply environmental error
        rho = errs.env_error_all(rho, 0, self.a1, time)

        N_ancillas = len(ancillas_pos)
        if N_ancillas == 1:
            N = len(rho.dims[0])
            p_success = ops.p_measurement_single_Zbasis(rho, projections[0],
                                                        N, ancillas_pos[0])
        # EPL success probability
        elif N_ancillas == 2 and projections == [1, 1]:
            p_success = ops.p_success_epl(rho, N=4, ancillas_pos=ancillas_pos)
        else:
            p_success = 1

        # Collapse the qubits in parrallel
        # Sort list to be able to reduce dimension and keep track of positions
        ancillas_pos = sorted(ancillas_pos)
        for i in range(N_ancillas):
            pos = ancillas_pos[i] - i
            rho = self._collapse_single(rho, pos,
                                        projections[i], "Z")
        return p_success, rho

    def _apply_two_qubit_gates(self, rho, controls, targets, sigma):
        # Apply one local two qubit Control type of gates in parallel
        # on each node.
        N = len(rho.dims[0])
        # NOTE: Number of gates is 1 because gates are applied in parallel
        self.check["two_qubit_gate"] += 1
        time = self.time_lookup["two_qubit_gate"]
        self.check["time"] += time
        self.check["time1"] += time
        # Apply environmental error
        rho = errs.env_error_all(rho, 0, self.a1, time)

        gates = self._get_two_qubit_gates(N, controls, targets, sigma)
        for i in range(len(gates)):
            rho = errs.two_qubit_gate(rho, gates[i], self.pg, N, controls[i],
                                      targets[i])
        return rho

    def start_bell_pair(self, rho=None):
        """Start with a raw Bell pair using the single click protocol."""
        self._reset_check()
        time, bell = self._generate_bell_single_click()
        self.check["time"] += time
        self.check["time0"] += time
        return 1, self.check, bell

    def start_BK(self, rho=None):
        """Circuit that generates a Bell pair from the Barret-Kok protocol."""
        self._reset_check()
        time, bell = self._generate_bell_pair_BK()
        self.check["time"] += time
        self.check["time0"] += time
        return 1, self.check, bell

    def start_epl(self, rho=None):
        """Circuit that generates a Bell pair using the EPL protocol."""
        self._reset_check()
        # Generate two raw Bell pairs
        time1, bell1 = self._generate_bell_single_click()
        time2, bell2 = self._generate_bell_single_click()
        self.check["time"] += time1 + time2
        self.check["time0"] += time1 + time2
        # Apply cutoff here

        # Dephase pair 1
        rho = errs.env_error_all(bell1, self.a0, self.a1, time2)
        rho = self._swap_pair(rho, [0, 1])

        # Join state
        rho = qt.tensor(rho, bell2)

        # NOTE: Decide which state collapses.
        # Apply two qubit gates
        controls = [0, 1]
        targets = [2, 3]
        rho = self._apply_two_qubit_gates(rho, controls,
                                          targets, "X")

        # Measure ancillas in Z basis
        projections = [1] * 2
        p_success, rho = self._collapse_ancillas_Z(rho, targets, projections)
        return p_success, self.check, rho

    def add_bell_pair(self, rho):
        """Append a raw Bell pair to the state."""
        self._reset_check()
        time, bell = self._generate_bell_pair()
        # Apply environmental error
        rho = errs.env_error_all(rho, self.a0, self.a1, time)
        self.check["time"] += time
        self.check["time1"] += time
        # Bell state is attached at the end of the complete state
        rho = qt.tensor(rho, bell)
        return 1, self.check, rho

    def swap_pair(self, rho, pair):
        """Does not really make a SWAP, only applies the noise."""
        self._reset_check()
        # Apply the noise
        self._swap_pair(rho, pair)
        return 1, self.check, rho

    def two_qubit_gates(self, rho, controls, targets, sigma):
        """
        Apply one local two qubit Control type of gates in parallel
        on each node.
        """
        self._reset_check()

        rho = self._apply_two_qubit_gates(rho, controls, targets, sigma)

        return 1, self.check, rho

    def collapse_ancillas_X(self, rho, ancillas_pos, projections):
        """
        Measure the ancillas in the X basis in parallel in each node.
        Ancillas position need to be the last part of the state.
        """
        # Reset check
        self._reset_check()
        p_success, rho = self._collapse_ancillas_X(rho, ancillas_pos,
                                                   projections)
        return p_success, self.check, rho

    def collapse_ancillas_Z(self, rho, ancillas_pos, projections):
        """
        Measure the ancillas in the X basis in parallel in each node.
        Ancillas position need to be the last part of the state.
        """
        # Reset check
        self._reset_check()
        p_success, rho = self._collapse_ancillas_Z(rho, ancillas_pos,
                                                   projections)
        return p_success, self.check, rho

    def single_selection(self, rho, operation_qubits, sigma):
        """
        Single selection round.
        Uses 2 ancilla qubits.
        """
        # Reset number of steps counter
        self._reset_check()

        # Generate raw bell pair
        rho = self._append_epl(rho)
        N = len(rho.dims[0])

        # Apply two qubit gates
        controls = [N-1, N-2]
        self._swap_pair(rho, controls)
        rho = self._apply_two_qubit_gates(rho, controls,
                                          operation_qubits, sigma)

        # Measure ancillas in X basis
        projections = [0] * 2
        p_success, rho = self._collapse_ancillas_X(rho, controls, projections)
        return p_success, self.check, rho

    def double_selection(self, rho, operation_qubits, sigma):
        """
        Double selection round.
        Uses 4 ancilla qubits.
        """
        # Reset number of steps counter
        self._reset_check()

        # Generate two bell pairs
        rho = self._append_epl(rho)
        rho = self._append_epl(rho)
        N = len(rho.dims[0])

        # Apply first two qubit gates
        controls1 = [N-3, N-4]
        self._swap_pair(rho, controls1)
        rho = self._apply_two_qubit_gates(rho, controls1,
                                          operation_qubits, sigma)

        # Apply second set of gates
        controls2 = [N-1, N-2]
        # Swap noise
        self._swap_pair(rho, controls2)
        rho = self._apply_two_qubit_gates(rho, controls2, controls1, "Z")

        # Measure ancillas in X basis
        projections = [0] * 2
        p_success2, rho = self._collapse_ancillas_X(rho, controls2,
                                                    projections)
        # Swap noise
        self._swap_pair(rho, controls1)
        p_success1, rho = self._collapse_ancillas_X(rho, controls1,
                                                    projections)
        p_success = p_success1 * p_success2
        return p_success, self.check, rho
