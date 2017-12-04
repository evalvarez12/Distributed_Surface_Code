"""
Individual circuits involved in the different steps of purification protocols.

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

    def __init__(self, ps, pm, pg, pn, p_env):
        """Init function."""
        # Set the parameters to all faulty opearations
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.pn = pn
        self.p_env = p_env

        # Detector efficiency in entanglement generation
        self.eta = 0.8

        # Dictionary to carry the check of all the operations made so far
        self.check = collections.Counter({"bell_pair": 0,
                                          "two_qubit_gate": 0,
                                          "single_qubit_gate": 0,
                                          "measurement": 0,
                                          "time": 0})

        # Lookup tables for the time it takes to make each operaation
        self.time_lookup = {"bell_pair": 1,
                            "two_qubit_gate": 1,
                            "single_qubit_gate": 1,
                            "measurement": 1}

    def change_parameters(self, ps, pm, pg, pn, p_env):
        """Function to change the parameters of all the operations."""
        # Reset all the parameters for the faulty operations
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.pn = pn
        self.p_env = p_env

        self._reset_check()

    def _reset_check(self):
        # Reset check dictionary
        self.check = collections.Counter(dict.fromkeys(self.check, 0))

    def _generate_bell_pair(self):
        # Probaility of success
        p_success = self.eta*(4 - self.eta)/4

        # This circuit number of steps and time it took
        attempts = self._success_number_of_attempts(p_success) + 1
        time = self.time_lookup["bell_pair"] * attempts

        # Update check
        self.check["bell_pair"] += 1
        # self.check["time"] += time

        # Generate noisy bell pair
        bell = errs.bell_pair(self.pn)
        return time, bell

    def _generate_bell_pair_BK(self):
        # Probaility of success
        p_success = self.eta**2/2

        # This circuit number of steps
        # Factor of 2 because uses twice the number of operations
        attempts = (self._success_number_of_attempts(p_success) + 1) * 2
        time = self.time_lookup["bell_pair"] * attempts

        # Update check
        self.check["bell_pair"] += 1

        # Generate noisy bell pair
        bell = errs.bell_pair(0.3*self.pn)
        return time, bell

    def _generate_noisy_plus(self):
        # This circuit time
        # NOTE:
        time = 1

        # Generate noisy plus
        plus = qt.snot() * qt.basis(2, 0)
        plus = plus * plus.dag()
        minus = qt.snot() * qt.basis(2, 1)
        minus = minus * minus.dag()
        plus = (1 - self.pm) * plus + self.pm * minus

        return time, plus

    def _success_number_of_attempts(self, p_success):
        # Up to 20 tries for success
        i = np.arange(100)
        d = self._distribution(p_success, i)
        return np.random.choice(i, 1, p=d)[0]
        # return 0

    def _distribution(self, p, n):
        # Distribution for the probability in the number of tries of
        # each event
        return p*(1-p)**n

    def _append_noisy_plus(self, rho):
        # Generate noisy plus
        time, plus = self.generate_noisy_plus()

        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, time, True)
        # Noisy plus tate is attached at the end of the complete state
        rho = qt.tensor(rho, plus)
        return rho

    def _append_bell_pair(self, rho):
        # Generate raw Bell pair to the state.
        time, bell = self._generate_bell_pair()

        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, time, True)

        # Bell state is attached at the end of the complete state
        rho = qt.tensor(rho, bell)
        return rho

    def _get_two_qubit_gates(self, N, controls, targets, sigma):
        # Construct two qubit Control type of gates.
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
        N = len(rho.dims[0])
        # NOTE: use only two CNOTs to perform a SWAP
        # Swap noise is only single qubit gate because one of the states
        # because a one way Swap gate is used

        # Apply noise
        for i in range(2):
            rho = errs.single_qubit_gate_noise(rho, self.pg, N, pos)
        return rho

    def _swap_pair(self, rho, pair):
        # NOTE: use only two CNOTs to perform a SWAP
        self.check["two_qubit_gate"] += 2
        self.check["time"] += self.time_lookup["two_qubit_gate"]*2
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
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, time, True)

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

    def _apply_two_qubit_gates(self, rho, controls, targets, sigma):
        # Apply one local two qubit Control type of gates in parallel
        # on each node.
        N = len(rho.dims[0])
        # NOTE: Number of gates is 1 because gates are applied in parallel
        self.check["two_qubit_gate"] += 1
        time = self.time_lookup["two_qubit_gate"]
        self.check["time"] += time

        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, time, True)

        gates = self._get_two_qubit_gates(N, controls, targets, sigma)
        for i in range(len(gates)):
            rho = errs.two_qubit_gate(rho, gates[i], self.pg, N, controls[i],
                                      targets[i])
        return rho

    def start_bell_pair(self, rho=None):
        """Start with a raw Bell pair to the state."""
        self._reset_check()
        time, bell = self._generate_bell_pair()
        self.check["time"] += time
        return 1, self.check, bell

    def add_bell_pair(self, rho):
        """Append a raw Bell pair to the state."""
        self._reset_check()
        time, bell = self._generate_bell_pair()
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, time, True)
        self.check["time"] += time
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

        # Swap noise
        self._swap_pair(rho, controls)

        rho = self._apply_two_qubit_gates(rho, controls, targets, sigma)

        return 1, self.check, rho

    def collapse_ancillas(self, rho, ancillas_pos, projections):
        """
        Measure the ancillas in the X basis in parallel in each node.
        Ancillas position need to be the last part of the state.
        """
        # Reset check
        self._reset_check()
        p_success, rho = self._collapse_ancillas_X(rho, ancillas_pos,
                                                   projections)
        return p_success, self.check, rho

    def single_selection(self, rho, operation_qubits, sigma):
        """
        Single selection round.
        Uses 2 ancilla qubits.
        """
        # Reset number of steps counter
        self._reset_check()

        # Swap noise
        self._swap_pair(rho, operation_qubits)

        # Generate raw bell pair
        rho = self._append_bell_pair(rho)
        N = len(rho.dims[0])

        # Apply two qubit gates
        controls = [N-1, N-2]
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

        # Swap noise
        self._swap_pair(rho, operation_qubits)

        # Generate two bell pairs
        rho = self._append_bell_pair(rho)
        rho = self._append_bell_pair(rho)
        N = len(rho.dims[0])

        # Apply first two qubit gates
        controls1 = [N-3, N-4]
        rho = self._apply_two_qubit_gates(rho, controls1,
                                          operation_qubits, sigma)

        # Apply second set of gates
        controls2 = [N-1, N-2]
        # Swap noise
        self._swap_pair(rho, controls1)
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
        print(p_success1, p_success2)
        return p_success, self.check, rho
