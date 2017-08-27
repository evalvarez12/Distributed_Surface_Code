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
            rho = errs.two_qubit_gate(rho, gates[i], self.pg, N, controls[i], targets[i])
        return rho

    def measure_single(self, rho, N, pos, basis):
        """
        Measure a single qubit in the state.
        Dimension is reduced after collapse
        """
        if basis == "X":
            measurement, collapsed_state = errs.measure_single_Xbasis(rho, self.pm, N, pos)
        elif basis == "Z":
            measurement, collapsed_state = errs.measure_single_Zbasis(rho, self.pm, N, pos)
        return measurement, collapsed_state


    # def measurements_Xbasis(self, rho, N, positions):
    #     measurements = []
    #     system_size = N
    #     for pos in positions:
    #         m, rho = self.measure_single(rho, system_size, pos, "X")
    #         system_size -= 1
    #         measurements += [m]
    #     return measurements, rho

    def operational_state(self, N):
        state = qt.bell_state('00')
        for i in range(1, N):
            state = qt.tensor(qt.bell_state('00'), state)
        state =  1/(2**N) * state * state.dag()
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
        rho = self.apply_two_qubit_gates(rho, N, controls, operation_qubits, sigma)

        # Measure ancillas in X basis
        success, collapsed_rho = self.collapse_check_success(rho, N, N_ancillas)
        return success, collapsed_rho

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
        rho = self.apply_two_qubit_gates(rho, N, controls, operation_qubits, sigma)

        # Apply second set of gates
        controls = [N-1, N-2]
        targets = [N-3, N-4]
        rho = self.apply_two_qubit_gates(rho, N, controls, targets, "Z")

        # Measure ancillas in X basis
        success, collapsed_rho = self.collapse_check_success(rho, N, N_ancillas)
        return success, collapsed_rho

    def one_dot(self, rho_initial, operation_qubits, sigma):
        """
        Perform the one dot procedure.
        Uses 4 ancillas.
        """
        N = len(rho_initial.dims[0]) + 2
        N_ancillas = 2

        success = False
        while not success:
            # Generate a raw Bell pair
            rho = self.append_bell_pair(rho_initial)

            # Rounds of single selection
            success, rho = self.single_selection(rho, [N-1, N-2], "X")
            if not success:
                continue
            success, rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Apply CNOT gates
        controls = [N-1, N-2]
        rho = self.apply_two_qubit_gates(rho, N, controls, operation_qubits, sigma)

        # Measure this procedures ancillas
        success, collapsed_rho = self.collapse_check_success(rho, N, N_ancillas)
        return success, collapsed_rho


    def two_dots(self, rho_initial, operation_qubits, sigma):
        """
        Perform the two dots procedure.
        Uses 4 ancillas.
        """
        N = len(rho_initial.dims[0]) + 2
        N_ancillas = 2

        success = False
        while not success:
            # Generate a raw Bell pair
            rho = self.append_bell_pair(rho_initial)

            # Rounds of single selection
            success, rho = self.single_selection(rho, [N-1, N-2], "X")
            if not success:
                continue
            success, rho = self.single_selection(rho, [N-1, N-2], "Z")

        # Apply CNOT gates
        controls = [N-1, N-2]
        rho = self.apply_two_qubit_gates(rho, N, controls, operation_qubits, sigma)

        # Extra round of single selection
        success, rho = self.single_selection(rho, [N-1, N-2], "Z")
        if not success:
            return False, None

        # Measure this procedures ancillas
        success, collapsed_rho = self.collapse_check_success(rho, N, N_ancillas)
        return success, collapsed_rho


    def expedient(self, stabilizer):
        """
        Perform the expedient protocol.
        Uses 4 data qubits and 12 ancillas.
        """

        # Initial state
        rho_initial = self.operational_state(4)
        rho_initial = rho_initial * rho_initial.dag()

        success = False
        while not success:
            print("--------->RESET")
            # Phase 1
            # First pair Bell state purification
            rho = self.append_bell_pair(rho_initial)
            N = len(rho.dims[0])
            operational_ancillas = [N-1, N-2]
            success, rho = self.double_selection(rho, operational_ancillas, "Z")
            if not success:
                continue
            success, rho = self.double_selection(rho, operational_ancillas, "X")
            if not success:
                continue

            # Second pair Bell state purification
            rho = self.append_bell_pair(rho)
            N = len(rho.dims[0])
            operational_ancillas = [N-1, N-2]
            success, rho = self.double_selection(rho, operational_ancillas, "Z")
            if not success:
                continue
            success, rho = self.double_selection(rho, operational_ancillas, "X")
            if not success:
                continue

            # Phase 2
            # Define pairs to form the GHZ state
            pair1 = [N-1, N-3]
            pair2 = [N-2, N-4]
            # Pair 1 operations
            success, rho = self.one_dot(rho, pair1, "Z")
            if not success:
                continue
            success, rho = self.one_dot(rho, pair1, "Z")
            if not success:
                continue

            #Pair 2 operations
            success, rho = self.one_dot(rho, pair2, "Z")
            if not success:
                continue
            success, rho = self.one_dot(rho, pair2, "Z")
            if not success:
                continue


        # Phase 3
        # Apply two qubit gates
        controls = [N-1, N-2, N-3, N-4]
        targets = [0, 1, 2, 3]
        rho = self.apply_two_qubit_gates(rho, N, controls, targets, stabilizer)
        measurements, rho = self.collapse_ancillas(rho, N, N_ancillas=4)
        return measurements, rho

    def stringent(self, stabilizer):
        """
        Perform the stringent protocol.
        Uses 4 data qubits and 12 ancillas.
        """

        # Initial state
        rho_initial = self.operational_state(4)
        rho_initial = rho_initial * rho_initial.dag()

        success = False
        while not success:
            print("--------->RESET")
            # Phase 1
            # First pair Bell state purification
            rho = self.append_bell_pair(rho_initial)
            N = len(rho.dims[0])
            operational_ancillas = [N-1, N-2]
            success, rho = self.double_selection(rho, operational_ancillas, "Z")
            if not success:
                continue
            success, rho = self.double_selection(rho, operational_ancillas, "X")
            if not success:
                continue
            success, rho = self.two_dot(rho, operational_ancillas, "Z")
            if not success:
                continue
            success, rho = self.two_dot(rho, operational_ancillas, "X")
            if not success:
                continue

            # Second pair Bell state purification
            rho = self.append_bell_pair(rho)
            N = len(rho.dims[0])
            operational_ancillas = [N-1, N-2]
            success, rho = self.double_selection(rho, operational_ancillas, "Z")
            if not success:
                continue
            success, rho = self.double_selection(rho, operational_ancillas, "X")
            if not success:
                continue
            success, rho = self.two_dot(rho, operational_ancillas, "Z")
            if not success:
                continue
            success, rho = self.two_dot(rho, operational_ancillas, "X")
            if not success:
                continue

            # Phase 2
            # Define pairs to form the GHZ state
            pair1 = [N-1, N-3]
            pair2 = [N-2, N-4]
            # Pair 1 operations
            success, rho = self.two_dot(rho, pair1, "Z")
            if not success:
                continue
            success, rho = self.two_dot(rho, pair1, "Z")
            if not success:
                continue

            #Pair 2 operations
            success, rho = self.two_dot(rho, pair2, "Z")
            if not success:
                continue
            success, rho = self.two_dot(rho, pair2, "Z")
            if not success:
                continue


        # Phase 3
        # Apply two qubit gates
        controls = [N-1, N-2, N-3, N-4]
        targets = [0, 1, 2, 3]
        rho = self.apply_two_qubit_gates(rho, N, controls, targets, stabilizer)
        measurements, rho = self.collapse_ancillas(rho, N, N_ancillas=4)
        return measurements, rho
