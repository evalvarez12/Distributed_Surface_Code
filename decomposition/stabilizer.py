"""
Routines for measuring stabilizers in a local system or in a disitributed
one by means of GHZ states.

author: Eduardo Villasenor
created-on: 05/02/18
"""
import qutip as qt
import numpy as np
import itertools
import error_models as errs
import tools.projectors as proj


class Stabilizer:
    """
    Class for holding functions relevant to making stabilizer measurements
    """

    def __init__(self, ps, pm, pg):
        """
        Init function.

        Parameters
        ----------
        ps : (scalar) single qubit gate error rate.
        pm : (scalar) measurement error rate.
        pg : (scalar) two qubit gate error rate.
        """
        self.ps = ps
        self.pm = pm
        self.pg = pg

    def change_parameters(self, ps, pm, pg):
        """
        Function to change the parameters of all the operations.

        Parameters
        ----------
        ps : (scalar) single qubit gate error rate.
        pm : (scalar) measurement error rate.
        pg : (scalar) two qubit gate error rate.
        """
        self.ps = ps
        self.pm = pm
        self.pg = pg

    def generate_noisy_plus(self):
        """Generate single noisy qubit in the |+> state."""
        plus = qt.snot() * qt.basis(2, 0)
        plus = plus * plus.dag()
        minus = qt.snot() * qt.basis(2, 1)
        minus = minus * minus.dag()
        plus = (1 - self.pm) * plus + self.pm * minus
        return plus

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

        Parameters
        ----------
        rho : (densmat) density matrix
        N : (int) total system size
        controls : (list) list with the control qubits position
        targets : (list) list with the target qubits position
        sigma : (string) X or Z depending on the gate
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

        Parameters
        ----------
        rho : (densmat) density matrix
        N : (int) total system size
        N_ancillas : (int) number of ancillas to be measured
        projections : (list) list with the projections in which every ancilla collapses
        """
        # The secuencial probabilities of finding the forced state
        probabilities = []
        for i in range(N_ancillas):
            rho = self.measure_single_forced(rho, N - i, N - i - 1,
                                             projections[i], "X")
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

    def measure_ghz_stabilizer_3on4(self, rho, ghz, parity_targets, stabilizer):
        """
        Measure a stabilizer through a GHZ state for the case when the size
        of the GHZ is one less than the number of qubits measured.

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
        if N_ghz != 3:
            raise ValueError("Measure stabilizer 3on4 dimension error")

        # Controls are the last qubits in rho
        controls = list(range(N - N_ghz, N)) + [N - 1]
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

    def local_stabilizer(self, rho, parity_targets, stabilizer):
        """
        Perform a local stabilizer measurement.
        Uses 4 data qubits and 1 ancillas.

        Parameters
        -----------
        rho : (densmat) state to be measured by the stabilizer
        parity_targets : (list) positions of the qubits to be measured
        stabilizer : (string) star or plaquette
        """
        # Append the initialized ancilla to the state|
        # NOTE: Initialized ancilla in mixed state
        ancilla = self.generate_noisy_plus()
        rho = qt.tensor(rho, ancilla)
        N = len(rho.dims[0])
        N_ancillas = 1

        # Apply two qubit gates
        controls = [N-1]*len(parity_targets)
        rho = self.apply_two_qubit_gates(rho, N, controls,
                                         parity_targets, stabilizer)

        # Get both projections of the state
        rho_even = self.collapse_ancillas_forced(rho, N, N_ancillas, [0])
        rho_odd = self.collapse_ancillas_forced(rho, N, N_ancillas, [1])

        # Get probabilites for each outcome
        p_even, p_odd = self._get_probabilities_measurement(rho, 1)
        return [p_even, p_odd], [rho_even, rho_odd]
