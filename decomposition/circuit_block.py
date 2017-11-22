"""
Individual circuits involved in the different steps of purification protocols.

author: Eduardo Villasenor
created-on: 20/11/17
"""
import qutip as qt
import numpy as np
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
        # Set the parameters to all faulty opearations
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.pn = pn
        self.p_env = p_env
        # Set number of steps of a block to 0
        self.n_steps = 0

    def change_parameters(self, ps, pm, pg, pn, p_env):
        """Function to change the parameters of all the operations."""
        # Reset all the parameters for the faulty operations
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.pn = pn
        self.p_env = p_env
        # Reset the number of steps
        self.n_steps = 0

    def generate_bell_pair(self):
        # This circuit number of steps
        steps = 1

        # Generate noisy bell pair
        bell = errs.bell_pair(self.pn)
        return steps, bell

    def generate_noisy_plus(self):
        # This circuit number of steps
        steps = 1

        # Generate noisy plus
        plus = qt.snot() * qt.basis(2, 0)
        plus = plus * plus.dag()
        minus = qt.snot() * qt.basis(2, 1)
        minus = minus * minus.dag()
        plus = (1 - self.pm) * plus + self.pm * minus

        return steps, plus

    def _append_noisy_plus(self, rho):
        plus = steps, self.generate_noisy_plus()

        # This circuit number of steps
        self.n_steps += steps
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, steps, True)
        # Noisy plus tate is attached at the end of the complete state
        rho = qt.tensor(rho, plus)
        return rho

    def _append_bell_pair(self, rho):
        # Append a raw Bell pair to the state.
        steps, bell = self.generate_bell_pair()

        # This circuit number of steps
        self.n_steps += steps
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, steps, True)

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
            collapsed_state = errs.measure_single_Xbasis_forced(rho,
                                                                self.pm,
                                                                project,
                                                                N,
                                                                pos)
        elif basis == "Z":
            collapsed_state = errs.measure_single_Zbasis_forced(rho,
                                                                self.pm,
                                                                project,
                                                                N,
                                                                pos)
        return collapsed_state

    def _collapse_ancillas_X(self, rho, ancillas_pos, projections):
        # Measure the ancillas in the X basis in parallel in each node.
        # NOTE: Ancillas position need to be the last part of the state.

        # Number of steps is divided over two because two measurements are made
        # in parallel over each node
        N_ancillas = len(ancillas_pos)
        steps = int(N_ancillas/2)
        self.n_steps += steps
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, steps, True)

        # Collapse the qubits in parrallel
        # Sort list to be able to reduce dimension and keep track of positions
        ancillas_pos = sorted(ancillas_pos)
        for i in range(N_ancillas):
            pos = ancillas_pos[i] - i
            rho = self._collapse_single(rho, pos,
                                        projections[i], "X")

        return rho

    def _apply_two_qubit_gates(self, rho, N, controls, targets, sigma):
        # Apply one local two qubit Control type of gates in parallel
        # on each node.

        # NOTE: This circuit number of steps is 1 because gates are applied in parallel
        steps = 1
        self.n_steps += steps
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, steps, True)

        gates = self._get_two_qubit_gates(N, controls, targets, sigma)
        for i in range(len(gates)):
            rho = errs.two_qubit_gate(rho, gates[i], self.pg, N, controls[i],
                                      targets[i])
        return rho

    def add_bell_pair(self, rho):
        """Append a raw Bell pair to the state."""
        steps, bell = self.generate_bell_pair()
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, steps, True)

        # Bell state is attached at the end of the complete state
        rho = qt.tensor(rho, bell)
        return 1, steps, rho

    def two_qubit_gates(self, rho, controls, targets, sigma):
        """
        Apply one local two qubit Control type of gates in parallel
        on each node.
        """
        N = len(rho.dims[0])
        # NOTE: This circuit number of steps is 1 because gates are applied in parallel
        steps = 1

        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, steps, True)

        gates = self._get_two_qubit_gates(N, controls, targets, sigma)
        for i in range(len(gates)):
            rho = errs.two_qubit_gate(rho, gates[i], self.pg, N, controls[i],
                                      targets[i])
        return 1, steps, rho

    def collapse_ancillas(self, rho, ancillas_pos, projections):
        """
        Measure the ancillas in the X basis in parallel in each node.
        Ancillas position need to be the last part of the state.
        """
        # Number of steps is  because it is made in parrallel
        # This circuit number of steps
        N_ancillas = len(ancillas_pos)
        steps = int(N_ancillas/2)
        N = len(rho.dims[0])
        # Apply environmental error
        rho = errs.env_dephasing_all(rho, self.p_env, steps, True)

        # Calculate probability of success using the same considerations as
        # in single selection
        p_success = ops.p_success_single_sel(rho, N, ancillas_pos)

        # Collapse the qubits in parrallel
        # Sort list to be able to reduce dimension and keep track of positions
        ancillas_pos = sorted(ancillas_pos)
        for i in range(N_ancillas):
            pos = ancillas_pos[i] - i
            rho = self._collapse_single(rho, pos,
                                        projections[i], "X")

        return p_success, steps, rho

    def single_selection(self, rho, operation_qubits, sigma):
        """
        Single selection round.
        Uses 2 ancilla qubits.
        """
        # Reset number of steps counter
        self.n_steps = 0
        # Generate raw bell pair
        rho = self._append_bell_pair(rho)
        N = len(rho.dims[0])
        N_ancillas = 2

        # Apply two qubit gates
        controls = [N-1, N-2]
        rho = self._apply_two_qubit_gates(rho, N, controls,
                                         operation_qubits, sigma)

        # Calculate probability of success
        p_success = ops.p_success_single_sel(rho, N, controls)

        # Measure ancillas in X basis
        projections = [0]*N_ancillas
        collapsed_rho = self._collapse_ancillas_X(rho, controls, projections)
        return p_success, self.n_steps, collapsed_rho

    def double_selection(self, rho, operation_qubits, sigma):
        """
        Double selection round.
        Uses 4 ancilla qubits.
        """
        # Reset number of steps counter
        self.n_steps = 0
        # Generate two bell pairs
        rho = self._append_bell_pair(rho)
        rho = self._append_bell_pair(rho)
        N = len(rho.dims[0])
        # Number of ancillas is 4, we asume all ancillas can be measured in parrallel
        N_ancillas = 4

        # Apply first two qubit gates
        controls1 = [N-3, N-4]
        rho = self._apply_two_qubit_gates(rho, N, controls1,
                                         operation_qubits, sigma)

        # Apply second set of gates
        controls2 = [N-1, N-2]
        rho = self._apply_two_qubit_gates(rho, N, controls2, controls1, "Z")

        # Calculate probability of success
        p_success = ops.p_success_double_sel(rho, N, controls1, controls2)

        # Measure ancillas in X basis
        ancillas_pos = list(range(N-1))
        projections = [0] * N_ancillas
        ancillas_pos = controls1 + controls2
        rho = self._collapse_ancillas_X(rho, ancillas_pos, projections)
        return p_success, self.n_steps, rho

    # def two_dots(self, rho, operation_qubits, sigma):
    #     """
    #     Hybrid purification round using single selection.
    #     Uses 4 ancilla qubits.
    #     """
    #     # Generate a bell pair and count the step
    #     rho = self._append_bell_pair(rho)
    #     steps = 1
    #     N = len(rho.dims[0])
    #     controls = [N-2, N-1]
    #
    #     # Rounds of single selection
    #     # Count steps manually because single sel resets self.steps
    #     p1, n1, rho = self.single_selection(rho, controls, "X")
    #     p2, n2, rho = self.single_selection(rho, controls, "X")
    #     steps += n1 + n2
    #     p_success = p1 * p2
    #
    #
    #     # Apply two qubit gates
    #     rho = self.apply_two_qubit_gates(rho, N, controls,
    #                                      operation_qubits, sigma)
    #     steps += 1
    #
    #     # Another round of single selection
    #     p, n, rho = self.single_selection(rho, controls, "Z")
    #     steps += n
    #     p_success *= p
    #
    #     # Calculate probability of success
    #     p = ops.p_success_single_sel(rho, N, controls)
    #     p_success *= p
    #
    #     # Measure ancillas in X basis
    #     projections = [0] * N_ancillas
    #     steps += 1
    #     rho = self.collapse_ancillas_forced(rho, N, N_ancillas, projections)
    #     return p_success, steps, rho
