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

    def __init__(self, ps, pm, pg, pn, system_size):
        self.ps = ps
        self.pm = pm
        self.pg = pg
        self.pn = pn
        self.N = system_size

    def generate_bell_pair(self):
        bell = errs.bell_pair(self.pn)
        bell = bell * bell.dag()
        return bell

    def append_bell_pair(self, rho):
        bell = self.generate_bell_pair()
        full_rho = qt.tensor(rho, bell)
        return full_rho

    def get_two_qubit_gates(self, N, controls, targets, sigma):
        gates = []
        if sigma == "X":
            for i in range(len(controls)):
                gates += [qt.cnot(N, controls[i], targets[i])]

        if sigma == "Z":
            for i in range(len(controls)):
                gates += [qt.cphase(np.pi, N, controls[i], targets[i])]

        return gates

    def apply_gates(self, rho, gates):
        for i in gates:
            rho = i * rho * i.dag()
        return rho

    def single_selection(self, rho, data_qubits, sigma):
        """
        Single selection round.
        Uses 4 qubits minimum.
        """
        # Generate raw bell pair
        rho = self.append_bell_pair(rho)
        N = len(rho.dims[0])

        # Apply two qubit gates
        controls = [N-1, N-2]
        gates = self.get_two_qubit_gates(N, controls, data_qubits, sigma)
        rho = self.apply_gates(rho, gates)

        # Measure ancillas in X basis
        measurement1, collapsed_state = ops.measure_single_Xbasis(rho, N, N-1, True)
        measurement2, collapsed_state = ops.measure_single_Xbasis(collapsed_state, N-1, N-2, True)
        return collapsed_state

    def double_selection(self, rho, data_qubits, sigma):
        """
        Double selection round.
        Uses 6 qubits minimum.
        """
        # Generate first bell pair
        rho = self.append_bell_pair(rho)
        N = len(rho.dims[0])

        # Apply two qubit gates
        controls = [N-1, N-2]
        gates = self.get_two_qubit_gates(N, controls, data_qubits, sigma)
        rho = self.apply_gates(rho, gates)

        # Generate second bell pair
        rho = self.append_bell_pair(rho)
        N = len(rho.dims[0])

        # Apply gates
        controls = [N-1, N-2]
        targets = [N-3, N-4]
        gates = self.get_two_qubit_gates(N, controls, data_qubits, "Z")
        rho = self.apply_gates(rho, gates)

        # Measure ancillas in X basis
        measurement1, collapsed_state = ops.measure_single_Xbasis(rho, N, N-1, True)
        measurement2, collapsed_state = ops.measure_single_Xbasis(collapsed_state, N-1, N-2, True)
        measurement3, collapsed_state = ops.measure_single_Xbasis(collapsed_state, N-2, N-3, True)
        measurement4, collapsed_state = ops.measure_single_Xbasis(collapsed_state, N-3, N-4, True)
        return collapsed_state
