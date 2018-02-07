"""
Routines for creating the protocols: EXPEDIENT, STRINGENT and MONOLITHIC
"""
import qutip as qt
import numpy as np
import itertools
import error_models as errs


class Tools:
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

    def change_parameters(self, ps, pm, pg, pn):
        """Function to change the parameters of all the operations."""
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.pn = pn

    def generate_bell_pair(self):
        """Generate a raw Bell pair."""
        bell = errs.bell_pair(self.pn)
        return bell

    def generate_noisy_plus(self):
        """Generate single noisy qubit in the |+> state."""
        plus = qt.snot() * qt.basis(2, 0)
        plus = plus * plus.dag()
        minus = qt.snot() * qt.basis(2, 1)
        minus = minus * minus.dag()
        plus = (1 - self.pm) * plus + self.pm * minus
        return plus

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

    def measure_single_forced(self, rho, N, pos, project, basis):
        """
        Measure a single qubit in the state.
        Dimension is reduced after collapse
        """
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

    def collapse_ancillas_forced(self, rho, N, N_ancillas, projections):
        """
        Measure the ancillas in the X basis.
        Ancillas position need to be the last part of the state.
        """
        # The secuencial probabilities of finding the forced state
        probabilities = []
        for i in range(N_ancillas):
            rho = self.measure_single_forced(rho, N - i, N - i - 1,
                                                projections[i], "X")


        return rho

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
        collapsed_rho = self.collapse_ancillas_forced(rho, N,
                                                             N_ancillas,
                                                             projections)
        return collapsed_rho

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
        collapsed_rho = self.collapse_ancillas_forced(rho, N,
                                                             N_ancillas,
                                                             projections)
        return collapsed_rho

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
        rho = self.single_selection(rho, [N-1, N-2], "X")
        rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Apply CNOT gates
        controls = [N-1, N-2]
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         operation_qubits,
                                         sigma)

        # Measure this procedures ancillas
        projections = [0] * N_ancillas
        collapsed_rho = self.collapse_ancillas_forced(rho, N,
                                                             N_ancillas,
                                                             projections)
        return collapsed_rho

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
        rho = self.single_selection(rho, [N-1, N-2], "X")
        rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Apply CNOT gates
        controls = [N-1, N-2]
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         operation_qubits,
                                         sigma)

        # Extra round of single selection
        rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Measure this procedures ancillas
        projections = [0] * N_ancillas
        collapsed_rho = self.collapse_ancillas_forced(rho, N,
                                                             N_ancillas,
                                                             projections)
        return collapsed_rho

    def make_ghz_expedient(self, N_ghz):
        """
        Perform the expedient protocol to generate a GHZ sate of four qubits.
        Uses 12 ancillas.
        """
        if N_ghz < 2:
            raise ValueError("Number of qubits in GHZ invalid")

        # Calculate number of bell pairs required
        N_pairs = np.int(N_ghz / 2)

        # Phase 1
        # Make first pair Bell state purification
        rho = self.generate_bell_pair()
        N = len(rho.dims[0])
        operational_ancillas = [N-1, N-2]
        rho = self.double_selection(rho, operational_ancillas, "Z")
        rho = self.double_selection(rho, operational_ancillas, "X")

        # Additional Bell states purification
        for i in range(N_pairs - 1):
            rho = self.append_bell_pair(rho)
            N = len(rho.dims[0])
            operational_ancillas = [N-1, N-2]
            rho = self.double_selection(rho, operational_ancillas, "Z")
            rho = self.double_selection(rho, operational_ancillas, "X")

        # Append single qubit if number of qubits is not pair
        if N_ghz % 2:
            rho = qt.tensor(rho, self.generate_noisy_plus())

        # Define pairs to form the GHZ state
        pairs = [[i, i + 2] for i in range(N_ghz - 2)]
        # Perform the pair operations
        for p in pairs:
            rho = self.one_dot(rho, p, "Z")
            rho = self.one_dot(rho, p, "Z")

        return rho

    def make_ghz_stringent(self, N_ghz):
        """
        Perform the expedient protocol to generate a GHZ sate of four qubits.
        Uses 12 ancillas.
        """
        if N_ghz < 2:
            raise ValueError("Number of qubits in GHZ invalid")

        # Calculate number of bell pairs required
        N_pairs = np.int(N_ghz / 2)

        # Phase 1
        # Make first pair Bell state purification
        rho = self.generate_bell_pair()
        N = len(rho.dims[0])
        operational_ancillas = [N-1, N-2]
        rho = self.double_selection(rho, operational_ancillas, "Z")
        rho = self.double_selection(rho, operational_ancillas, "X")
        rho = self.two_dots(rho, operational_ancillas, "Z")
        rho = self.two_dots(rho, operational_ancillas, "X")

        # Additional Bell states purification
        for i in range(N_pairs - 1):
            rho = self.append_bell_pair(rho)
            N = len(rho.dims[0])
            operational_ancillas = [N-1, N-2]
            _, rho = self.double_selection(rho, operational_ancillas, "Z")
            _, rho = self.double_selection(rho, operational_ancillas, "X")
            _, rho = self.two_dots(rho, operational_ancillas, "Z")
            _, rho = self.two_dots(rho, operational_ancillas, "X")

        # Append single qubit if number of qubits is not pair
        if N_ghz % 2:
            rho = qt.tensor(rho, self.generate_noisy_plus())

        # Define pairs to form the GHZ state
        pairs = [[i, i + 2] for i in range(N_ghz - 2)]
        # Perform the pair operations
        for p in pairs:
            rho = self.two_dots(rho, p, "Z")
            rho = self.two_dots(rho, p, "Z")

        return rho

    def twirl_ghz(self, ghz):
        N = len(ghz.dims[0])
        twirled_state = ghz * 0
        permutations = itertools.permutations(range(N), N)
        for p in permutations:
            twirled_state += ghz.permute(p)
        return twirled_state/np.math.factorial(N)

    def measure_ghz_stabilizer(self, rho_initial, ghz, parity_targets, stabilizer):
        # Apply two qubit gates
        N_ghz = len(ghz.dims[0])
        rho = qt.tensor(rho_initial, ghz)
        N = len(rho.dims[0])
        # Controls are the last qubits in rho
        controls = list(range(N - N_ghz, N))
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         parity_targets, stabilizer)
        projections_even = [0] * N_ghz
        rho_even = self.collapse_ancillas_forced(rho, N,
                                                         N_ghz,
                                                         projections_even)
        projections_odd = [0] * N_ghz
        projections_odd[-1] = 1
        rho_odd = self.collapse_ancillas_forced(rho, N,
                                                       N_ghz,
                                                       projections_odd)
        # NOTE: Cheating by putting probs by hand
        p_odd = 0.5
        p_even = 0.5
        return [p_even, p_odd], [rho_even, rho_odd]

    def measure_ghz_stabilizer_3on4(self, rho_initial, ghz, parity_targets, stabilizer):
        # Apply two qubit gates
        N_ghz = len(ghz.dims[0])
        rho = qt.tensor(rho_initial, ghz)
        N = len(rho.dims[0])
        if N_ghz != 3:
            raise ValueError("Measure stabilizer 3on4 dimension error")

        # Controls are the last qubits in rho
        controls = list(range(N - N_ghz, N)) + [N - 1]
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         parity_targets, stabilizer)
        projections_even = [0] * N_ghz
        rho_even = self.collapse_ancillas_forced(rho, N,
                                                 N_ghz,
                                                 projections_even)
        projections_odd = [0] * N_ghz
        projections_odd[-1] = 1
        rho_odd = self.collapse_ancillas_forced(rho, N,
                                                N_ghz,
                                                projections_odd)
        # NOTE: Cheating by putting probs by hand
        p_odd = 0.5
        p_even = 0.5
        return [p_even, p_odd], [rho_even, rho_odd]


    def expedient(self, rho_initial, parity_targets, stabilizer):
        """
        Perform the expedient protocol.
        Uses 4 data qubits and 12 ancillas.
        """
        # GHZ number of qubits is same as the number of qubits
        # in the state to be parity measured
        N_parity = len(parity_targets)
        ghz = self.make_ghz_expedient(N_parity)
        ghz = self.twirl_ghz(ghz)
        return self.measure_ghz_stabilizer(rho_initial, ghz,
                                           parity_targets,
                                           stabilizer)

    def stringent(self, rho_initial, parity_targets, stabilizer):
        """
        Perform the stringent protocol.
        Uses 4 ancillas per data qubit at maximum.
        """
        # GHZ number of qubits is same as the number of qubits
        # in the state to be parity measured
        N_parity = len(parity_targets)
        ghz = self.make_ghz_stringent(N_parity)
        ghz = self.twirl_ghz(ghz)
        return self.measure_ghz_stabilizer(rho_initial, ghz,
                                           parity_targets,
                                           stabilizer)

    def local_stabilizer(self, rho_initial, parity_targets, stabilizer):
        """
        Perform the monolithic stabilizer protocol.
        Uses 4 data qubits and 1 ancillas.
        """
        # Append the initialized ancilla to the state|
        # NOTE: Naomi starts adding error to the initialized
        ancilla = self.generate_noisy_plus()
        rho = qt.tensor(rho_initial, ancilla)
        N = len(rho.dims[0])
        N_ancillas = 1

        # Apply two qubit gates
        controls = [N-1]*len(parity_targets)
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         parity_targets, stabilizer)

        # Get both projections of the state
        rho_even = self.collapse_ancillas_forced(rho,
                                                         N,
                                                         N_ancillas,
                                                         [0])
        rho_odd = self.collapse_ancillas_forced(rho,
                                                       N,
                                                       N_ancillas,
                                                       [1])
        # NOTE: Cheating by putting probs by hand
        p_odd = 0.5
        p_even = 0.5
        return [p_even, p_odd], [rho_even, rho_odd]
