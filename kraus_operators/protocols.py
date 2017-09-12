"""
Routines for creating the protocols: EXPEDIENT, STRINGENT and MONOLITHIC
"""
import qutip as qt
import numpy as np
import error_models as errs
import operations as ops

class Protocols:
    """
    Class for holding all protocols.

    Paramenters
    -----------
    ps - single qubit gate error
    pm - single qubit measurement error
    pg - two-qubit gate error
    pn - network error
    """

    def __init__(self, ps, pm, pg, pn):
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.pn = pn

    def fidelity(self, rhoA, stateB):
        """
        Fidelity for the special case when one of the states is a pure state.
        """
        return (stateB.dag() * rhoA * stateB).norm()

    def generate_bell_pair(self):
        """Generate a raw Bell pair."""
        bell = errs.bell_pair(self.pn)
        bell = bell * bell.dag()
        return bell

    def append_bell_pair(self, rho):
        """Append a raw Bell pair to the state."""
        bell = self.generate_bell_pair()
        # Bell state is attached at the end of the complete state
        full_rho = qt.tensor(rho, bell)
        return full_rho

    def get_two_qubit_gates(self, N, controls, targets, sigma):
        """
        Construct two qubit Control type of gates.
        """
        gates = []
        if sigma == "X":
            for i in range(len(controls)):
                gates += [qt.cnot(N, controls[i], targets[i])]

        if sigma == "Z":
            for i in range(len(controls)):
                gates += [qt.cphase(np.pi, N, controls[i], targets[i])]

        return gates

    def apply_two_qubit_gates(self, rho, N, controls, targets, sigma):
        """
        Apply a series of two qubit Control type of gates.
        """
        gates = self.get_two_qubit_gates(N, controls, targets, sigma)
        for i in range(len(gates)):
            rho = errs.two_qubit_gate(rho, gates[i], self.pg, N, controls[i],
                                      targets[i])
        return rho

    def parity_projection_ket(self, psi, targets, measurement, parity):
        plus = qt.snot() * qt.basis(2, 0)
        full_psi = qt.tensor(psi, plus)
        N = len(full_psi.dims[0])
        control = N-1
        if parity == "X":
            for t in targets:
                full_psi = qt.cnot(N, control, t) * full_psi
        if parity == "Z":
            for t in targets:
                full_psi = qt.cphase(np.pi, N, control, t) * full_psi

        collapsed_psi = ops.collapse_single_Xbasis_ket(full_psi, measurement,
                                                       N, N-1, True)
        return collapsed_psi

    def measure_single(self, rho, N, pos, basis):
        """
        Measure a single qubit in the state.
        Dimension is reduced after collapse
        """
        if basis == "X":
            measurement, collapsed_state = errs.measure_single_Xbasis(rho,
                                                                      self.pm,
                                                                      N,
                                                                      pos)
        elif basis == "Z":
            measurement, collapsed_state = errs.measure_single_Zbasis(rho,
                                                                      self.pm,
                                                                      N,
                                                                      pos)
        return measurement, collapsed_state

    def measure_single_forced(self, rho, N, pos, project, basis):
        """
        Measure a single qubit in the state.
        Dimension is reduced after collapse
        """
        if basis == "X":
            p, collapsed_state = errs.measure_single_Xbasis_forced(rho,
                                                                   self.pm,
                                                                   project,
                                                                   N,
                                                                   pos)
        elif basis == "Z":
            p, collapsed_state = errs.measure_single_Zbasis_forced(rho,
                                                                   self.pm,
                                                                   project,
                                                                   N,
                                                                   pos)
        return p, collapsed_state

    def operational_state_ket(self, N):
        """
        Create a bipartite state.
        Also known as operational state.
        """
        state = qt.bell_state('00')
        for i in range(1, N):
            state = qt.tensor(qt.bell_state('00'), state)
        return state

    def collapse_ancillas(self, rho, N, N_ancillas):
        """
        Measure the ancillas in the X basis.
        Ancillas position need to be the last part of the state.
        """
        measurements = []
        for i in range(N_ancillas):
            m, rho = self.measure_single(rho, N - i, N - i - 1, "X")
            measurements += [m]

        return measurements, rho

    def collapse_ancillas_forced(self, rho, N, N_ancillas, projections):
        """
        Measure the ancillas in the X basis.
        Ancillas position need to be the last part of the state.
        """
        # The secuencial probabilities of finding the forced state
        probabilities = []
        for i in range(N_ancillas):
            p, rho = self.measure_single_forced(rho, N - i, N - i - 1,
                                                projections[i], "X")
            probabilities += [p]

        return probabilities, rho

    def collapse_check_success(self, rho, N, N_ancillas):
        """
        Measure the ancillas in the X basis.
        Ancillas position need to be the last part of the state.
        """
        measurements, rho = self.collapse_ancillas(rho, N, N_ancillas)
        if len(set(measurements)) > 1:
            return False, None

        return True, rho

    def single_selection(self, rho, operation_qubits, sigma):
        """
        Single selection round.
        Uses 2 ancilla qubits.
        """
        # Generate raw bell pair
        rho = self.append_bell_pair(rho)
        N = len(rho.dims[0])
        N_ancillas = 2

        # Apply two qubit gates
        controls = [N-1, N-2]
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         operation_qubits, sigma)

        # Measure ancillas in X basis
        projections = [0]*N_ancillas
        probs, collapsed_rho = self.collapse_ancillas_forced(rho, N,
                                                             N_ancillas,
                                                             projections)
        return probs, collapsed_rho

    def double_selection(self, rho, operation_qubits, sigma):
        """
        Double selection round.
        Uses 4 ancilla qubits.
        """
        # Generate two bell pairs
        rho = self.append_bell_pair(rho)
        rho = self.append_bell_pair(rho)
        N = len(rho.dims[0])
        N_ancillas = 4

        # Apply first two qubit gates
        controls = [N-3, N-4]
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         operation_qubits, sigma)

        # Apply second set of gates
        controls = [N-1, N-2]
        targets = [N-3, N-4]
        rho = self.apply_two_qubit_gates(rho, N, controls, targets, "Z")

        # Measure ancillas in X basis
        projections = [0] * N_ancillas
        probs, collapsed_rho = self.collapse_ancillas_forced(rho, N,
                                                             N_ancillas,
                                                             projections)
        return probs, collapsed_rho

    def one_dot(self, rho_initial, operation_qubits, sigma):
        """
        Perform the one dot procedure.
        Uses 4 ancillas.
        """
        N = len(rho_initial.dims[0]) + 2
        N_ancillas = 2

        # Generate a raw Bell pair
        rho = self.append_bell_pair(rho_initial)

        # Rounds of single selection
        _, rho = self.single_selection(rho, [N-1, N-2], "X")
        _, rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Apply CNOT gates
        controls = [N-1, N-2]
        rho = self.apply_two_qubit_gates(rho, N, controls, operation_qubits,
                                         sigma)

        # Measure this procedures ancillas
        projections = [0] * N_ancillas
        probs, collapsed_rho = self.collapse_ancillas_forced(rho, N,
                                                             N_ancillas,
                                                             projections)
        return probs, collapsed_rho

    def two_dots(self, rho_initial, operation_qubits, sigma):
        """
        Perform the two dots procedure.
        Uses 4 ancillas.
        """
        N = len(rho_initial.dims[0]) + 2
        N_ancillas = 2

        # Generate a raw Bell pair
        rho = self.append_bell_pair(rho_initial)

        # Rounds of single selection
        _, rho = self.single_selection(rho, [N-1, N-2], "X")
        _, rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Apply CNOT gates
        controls = [N-1, N-2]
        rho = self.apply_two_qubit_gates(rho, N, controls, operation_qubits, sigma)

        # Extra round of single selection
        _, rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Measure this procedures ancillas
        projections = [0] * N_ancillas
        probs, collapsed_rho = self.collapse_ancillas_forced(rho, N, N_ancillas, projections)
        return probs, collapsed_rho

    def expedient(self, rho_initial, parity_targets, stabilizer):
        """
        Perform the expedient protocol.
        Uses 4 data qubits and 12 ancillas.
        """
        # Phase 1
        # First pair Bell state purification
        rho = self.append_bell_pair(rho_initial)
        N = len(rho.dims[0])
        operational_ancillas = [N-1, N-2]
        _, rho = self.double_selection(rho, operational_ancillas, "Z")
        _, rho = self.double_selection(rho, operational_ancillas, "X")

        # Second pair Bell state purification
        rho = self.append_bell_pair(rho)
        N = len(rho.dims[0])
        operational_ancillas = [N-1, N-2]
        _, rho = self.double_selection(rho, operational_ancillas, "Z")
        _, rho = self.double_selection(rho, operational_ancillas, "X")

        # Phase 2
        # Define pairs to form the GHZ state
        pair1 = [N-1, N-3]
        pair2 = [N-2, N-4]
        # Pair 1 operations
        _, rho = self.one_dot(rho, pair1, "Z")
        _, rho = self.one_dot(rho, pair1, "Z")

        #Pair 2 operations
        _, rho = self.one_dot(rho, pair2, "Z")
        _, rho = self.one_dot(rho, pair2, "Z")


        # Phase 3
        # Apply two qubit gates
        controls = [N-1, N-2, N-3, N-4]
        targets = parity_targets
        rho = self.apply_two_qubit_gates(rho, N, controls, targets, stabilizer)
        measurements, rho = self.collapse_ancillas(rho, N, N_ancillas=4)
        return measurements, rho

    def stringent(self, rho_initial, parity_targets, stabilizer):
        """
        Perform the stringent protocol.
        Uses 4 data qubits and 12 ancillas.
        """
        # Phase 1
        # First pair Bell state purification
        rho = self.append_bell_pair(rho_initial)
        N = len(rho.dims[0])
        operational_ancillas = [N-1, N-2]
        _, rho = self.double_selection(rho, operational_ancillas, "Z")
        _, rho = self.double_selection(rho, operational_ancillas, "X")
        _, rho = self.two_dots(rho, operational_ancillas, "Z")
        _, rho = self.two_dots(rho, operational_ancillas, "X")

        # Second pair Bell state purification
        rho = self.append_bell_pair(rho)
        N = len(rho.dims[0])
        operational_ancillas = [N-1, N-2]
        _, rho = self.double_selection(rho, operational_ancillas, "Z")
        _, rho = self.double_selection(rho, operational_ancillas, "X")
        _, rho = self.two_dots(rho, operational_ancillas, "Z")
        _, rho = self.two_dots(rho, operational_ancillas, "X")

        # Phase 2
        # Define pairs to form the GHZ state
        pair1 = [N-1, N-3]
        pair2 = [N-2, N-4]
        # Pair 1 operations
        _, rho = self.two_dots(rho, pair1, "Z")
        _, rho = self.two_dots(rho, pair1, "Z")

        # Pair 2 operations
        _, rho = self.two_dots(rho, pair2, "Z")
        _, rho = self.two_dots(rho, pair2, "Z")

        # Phase 3
        # Apply two qubit gates
        controls = [N-1, N-2, N-3, N-4]
        targets = parity_targets
        rho = self.apply_two_qubit_gates(rho, N, controls, targets, stabilizer)
        measurements, rho = self.collapse_ancillas(rho, N, N_ancillas=4)
        return measurements, rho

    def monolithic(self, rho_initial, parity_targets, stabilizer):
        """
        Perform the monolithic stabilizer protocol.
        Uses 4 data qubits and 1 ancillas.
        """
        # Append the initialized ancilla to the state|
        # NOTE: Naomi starts adding error to the initialized
        plus = qt.snot() * qt.basis(2, 0)
        plus = plus * plus.dag()
        minus = qt.snot() * qt.basis(2, 1)
        minus = minus * minus.dag()
        ancilla = (1 - self.pm)*plus + self.pm*minus

        # ancilla = qt.snot() * qt.basis(2, 0)
        # ancilla = ancilla * ancilla.dag()
        rho = qt.tensor(rho_initial, ancilla)
        N = len(rho.dims[0])
        N_ancillas = 1

        # Apply two qubit gates
        controls = [N-1]*len(parity_targets)
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         parity_targets, stabilizer)

        # Get both projections of the state
        p_even, rho_even = self.collapse_ancillas_forced(rho,
                                                         N,
                                                         N_ancillas,
                                                         [0])
        p_odd, rho_odd = self.collapse_ancillas_forced(rho,
                                                       N,
                                                       N_ancillas,
                                                       [1])
        return [p_even[0], p_odd[0]], [rho_even, rho_odd]
