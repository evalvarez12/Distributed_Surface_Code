"""
Routines for creating the protocols of Naomi Nickerson:
EXPEDIENT, STRINGENT and MONOLITHIC

NOTE: This is not used anymore

author: Eduardo Villasenor
created-on: 12/08/17
"""
import qutip as qt
import numpy as np
import itertools
import error_models as errs
import tools.projectors as proj


class Protocols:
    """
    Class for holding GHZ formation protocolsprotocols.

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

    def raw_state(self, pn, r):
        imperfect_state = errs.bell_pair_psi(pn)
        noise = qt.tensor(qt.basis(2, 1), qt.basis(2, 1))
        noise = noise * noise.dag()
        return (1 - r) * imperfect_state + r * noise

    def generate_epl(self):
        """
        Circuit that generates a Bell pair using the EPL protocol.
        """
        # Define extra paramters
        r = 1/4.
        p_drift = 1/100.

        # Generate states
        raw1 = self.raw_state(self.pn, r)
        raw2 = self.raw_state(self.pn, r)
        raw2 = errs.drift(raw2, p_drift)

        # Join state
        rho = qt.tensor(raw1, raw2)

        # Apply two qubit gates
        controls = [0, 1]
        targets = [2, 3]
        rho = self.apply_two_qubit_gates(rho, 4, controls,
                                         targets, "X")

        # Collapse ancillas according to EPL
        N = 4
        pos = 2
        for i in range(2):
            rho = self.measure_single_forced(rho, N, pos,
                                             1, "Z")
            N = N - 1

        X = qt.rx(np.pi, 2, 0)
        rho = errs.single_qubit_gate(rho, X, self.ps, 2, 0)

        return rho

    def twirl_ghz(self, ghz):
        """
        Twirl the GHZ state to disitribute error uniformly.
        Twirling is made trough averaging through permutation of all qubits in the GHZ.

        Parameters
        -----------
        ghz : (densmat) density matrix of a ghz state.
        """
        N = len(ghz.dims[0])
        twirled_state = ghz * 0
        # Get all permuations and average over them
        permutations = itertools.permutations(range(N), N)
        for p in permutations:
            twirled_state += ghz.permute(p)
        return twirled_state/np.math.factorial(N)

    def measure_ghz_stabilizer(self, rho, ghz, parity_targets, stabilizer):
        """
        Measure a stabilizer through a GHZ state.
        Size of the GHZ equals the number of qubits measured through the stabilizer.

        Parameters
        -----------
        rho : (densmat) state to be measured by the stabilizer
        ghz : (densmat) GHZ state to be used in the stabilizer measurement
        parity_targets : (list) positions of the qubits to be measured
        stabilizer : (string) star or plaquette
        """
        # Apply two qubit gates
        N_ghz = len(ghz.dims[0])
        rho = qt.tensor(rho, ghz)
        N = len(rho.dims[0])
        # Controls are the last qubits in rho
        controls = list(range(N - N_ghz, N))
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         parity_targets, stabilizer)
        projections_even = [0] * N_ghz
        rho_even = self.collapse_ancillas_forced(rho, N, N_ghz,
                                                 projections_even)
        projections_odd = [0] * N_ghz
        projections_odd[-1] = 1
        rho_odd = self.collapse_ancillas_forced(rho, N, N_ghz,
                                                projections_odd)

        # Get probabilites for each outcome
        p_even, p_odd = self._get_probabilities_measurement(rho, N_ghz)
        return [p_even, p_odd], [rho_even, rho_odd]

    def _get_probabilities_measurement(self, rho, N_ghz):
        # Compute the probabilites of getting a even or a odd measurement result
        # Size of the state that is not the GHZ
        N = len(rho.dims[0]) - N_ghz
        # Get the sum of ven and odd projectors
        P_even = proj.even_projectors(N_ghz)
        P_odd = proj.odd_projectors(N_ghz)

        # Tensor with corresponding identity
        # Remember GHZ state is at the end of the state
        if N != 0:
            P_even = qt.tensor(qt.qeye([2]*N), P_even)
            P_odd = qt.tensor(qt.qeye([2]*N), P_odd)

        # Compute probabilites
        p_even = (P_even * rho * P_even.dag()).tr()
        p_odd = (P_odd * rho * P_odd.dag()).tr()
        return [p_even, p_odd]

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
        print("N_GHZ: ", N_ghz)
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
            rho = self.double_selection(rho, operational_ancillas, "Z")
            rho = self.double_selection(rho, operational_ancillas, "X")
            rho = self.two_dots(rho, operational_ancillas, "Z")
            rho = self.two_dots(rho, operational_ancillas, "X")

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

    def make_ghz_basic(self):
        N_ghz=4
        # Calculate number of bell pairs required
        N_pairs = np.int(N_ghz)

        # Phase 1
        # Make first pair Bell state
        rho = self.generate_epl()

        # Additional Bell states purification
        for i in range(N_pairs - 1):
            rho = qt.tensor(rho, self.generate_epl())

        N = len(rho.dims[0])
        N_ancillas = 4

        controls = [1, 3, 0, 2]
        targets = [4, 5, 6, 7]
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         targets,
                                         "Z")

        # Measure this procedures ancillas
        projections = [0] * len(controls)
        rho = self.collapse_ancillas_forced(rho, N,
                                            N_ancillas,
                                            projections)

        return rho

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

    def basic(self, rho_initial, parity_targets, stabilizer):
        """
        Perform the basic protocol.
        Uses 4 ancillas per data qubit at maximum.
        """
        # GHZ number of qubits is same as the number of qubits
        # in the state to be parity measured
        N_parity = len(parity_targets)
        ghz = self.make_ghz_basic()

        ghz = self.ghz_purification(ghz)

        ghz = self.twirl_ghz(ghz)
        return self.measure_ghz_stabilizer(rho_initial, ghz,
                                           parity_targets,
                                           stabilizer)

    def ghz_purification(self, ghz):
        # Two GHZ states to make a parity measuement
        ghz = qt.tensor(ghz, ghz)

        # Apply second set of gates
        targets = [0, 1, 2, 3]
        controls = [4, 5, 6, 7]
        ghz = self.apply_two_qubit_gates(ghz, 8, controls, targets, "X")

        # Collapse one GHZ
        # Measure this procedures ancillas
        projections = [0] * 4
        ghz = self.collapse_ancillas_forced(ghz, 8, 4,
                                            projections)
        return ghz
