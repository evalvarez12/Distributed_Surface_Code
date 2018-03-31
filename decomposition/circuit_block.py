"""
Individual circuits involved the remote generation of entanglement
and in different steps of purification protocols. Must be used together
with 'circuit.py' to assemeble complete circuits.

Includes features to keep track of how many operations and time has passed
after each circuit.

author: Eduardo Villasenor
created-on: 20/11/17
"""
import qutip as qt
import numpy as np
import collections
import tools.operations as ops
import error_models as errs
import tools.p_success as ps
import tools.ghz_correction as gc


class Blocks:
    """
    Class for holding all protocols.
    Each circuit block returns the resulting state, number of steps used
    and the probability of success if applicable.
    """

    def __init__(self, ps, pm, pg, eta, a0, a1, theta):
        """Init function.

        Parameters
        ----------
        ps : (scalar) single qubit gate error rate.
        pm : (scalar) measurement error rate.
        pg : (scalar) two qubit gate error rate.
        eta : (scalar) detection efficiency.
        a0 : (scalar) extra environmental error when electron spin is being operated.
        a1 : (scalar) default environmental error.
        theta : (scalar) determines how the states are initialized when generating remote
                entanglement.
        """
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
        """
        Function to change the parameters of all the operations.

        Parameters
        ----------
        ps : (scalar) single qubit gate error rate.
        pm : (scalar) measurement error rate.
        pg : (scalar) two qubit gate error rate.
        eta : (scalar) detection efficiency.
        a0 : (scalar) extra environmental error when electron spin is being operated.
        a1 : (scalar) default environmental error.
        theta : (scalar) determines how the states are initialized when generating remote
                entanglement.
        """
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
        # Reset to 0 all vaulues from check
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
        # For this case set the initial state to be |+>
        # Probaility of success
        s = .5
        r = (1 - self.eta)*s/(1 - self.eta*s)
        p_success = (1 - r)*self.eta**2

        # This circuit number of steps
        attempts = self._success_number_of_attempts(p_success) + 1
        time = self.time_lookup["bell_pair"] * attempts

        # Update check
        self.check["bell_pair"] += 1

        # Generate Bell pair with F=1
        bell = qt.bell_state('10') * qt.bell_state('10').dag()
        # Apply noisy 1q X to transform state
        bell = errs.single_qubit_gate(bell, qt.rx(np.pi, 2, 0), self.ps, 2, 0)

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

        # NOTE: This is expensive, adjust manually as required
        # 100000 tries required for EPL
        # 1000000 for BK
        i = np.arange(100000)
        d = self._distribution(p_success, i)
        return np.random.choice(i, 1, p=d)[0]

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
            rho = errs.measure_single_Xbasis_forced(rho, self.pm,
                                                    project, N, pos)
        elif basis == "Z":
            rho = errs.measure_single_Zbasis_forced(rho, self.pm,
                                                    project, N, pos)
        return rho

    def _swap_noise(self, rho, pos):
        # Apply the noise induced by the one way SWAP gate
        # NOTE: use only two CNOTs to perform a SWAP
        # Swap noise is only single qubit gate because one of the states
        # because a one way Swap gate is used
        ancilla = qt.basis(2, 0) * qt.basis(2, 0).dag()
        ancilla = errs.env_error_all(ancilla, self.a0, self.a1,
                                     self.check["time0"])
        rho = qt.tensor(rho, ancilla)
        N = len(rho.dims[0])

        # Define ideal gates
        CNOT1 = qt.cnot(N, N-1, pos)
        CNOT2 = qt.cnot(N, pos, N-1)

        # Apply 3 CNOTS to get a swap
        rho = errs.two_qubit_gate(rho, CNOT1, self.pg, N, pos, N-1)
        rho = errs.two_qubit_gate(rho, CNOT2, self.pg, N, pos, N-1)
        rho = errs.two_qubit_gate(rho, CNOT1, self.pg, N, pos, N-1)

        # Measure the ancilla to reduce the dimension
        # m, rho =  ops.random_measure_single_Zbasis(rho, N, pos, True)
        # print("M:", m)
        rho = ops.collapse_single_Zbasis(rho, 0, N, pos, True)
        return rho

    def _swap_pair(self, rho, pair):
        # Apply the noise due to SWAP on two states
        # NOTE: use only two CNOTs to perform a SWAP
        self.check["two_qubit_gate"] += 3
        time = self.time_lookup["two_qubit_gate"] * 3
        self.check["time"] += time
        self.check["time1"] += time
        rho = self._swap_noise(rho, pair[0])
        rho = self._swap_noise(rho, pair[1])
        rho = errs.env_error_all(rho, 0, self.a1, time)
        return rho

    def _collapse_ancillas_X(self, rho, ancillas_pos, projections):
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
        N = len(rho.dims[0])

        # Apply error due to the needed rotation, use I instead of real H gate
        identity = qt.qeye([2]*N)
        rho = self._apply_single_qubit_gates(rho, [identity, identity],
                                             ancillas_pos)

        # Collapse the qubits in parrallel
        # Sort list to be able to reduce dimension and keep track of positions
        ancillas_pos = sorted(ancillas_pos)
        for i in range(N_ancillas):
            pos = ancillas_pos[i] - i
            rho = self._collapse_single(rho, pos,
                                        projections[i], "X")

        return 1, rho

    def _collapse_ancillas_EPL(self, rho, ancillas_pos):
        # For the special EPL case
        # Measure the ancillas in the Z basis in parallel in each node.
        # NOTE: All ancillas are collapsed in parallel, Hadamard operations
        # are used to measure on X basis
        self.check["measurement"] += 1

        time = self.time_lookup["measurement"]
        self.check["time"] += time
        self.check["time1"] += time
        # Apply environmental error
        rho = errs.env_error_all(rho, 0, self.a1, time)

        N_ancillas = len(ancillas_pos)
        if N_ancillas != 2:
            raise ValueError("EPL ancillas collapse: Bad number of ancillas")
        p_success = ps.epl(rho, N=4, ancillas_pos=ancillas_pos)
        projections = [1, 1]

        # Collapse the qubits in parrallel
        # Sort list to be able to reduce dimension and keep track of positions
        ancillas_pos = sorted(ancillas_pos)
        for i in range(N_ancillas):
            pos = ancillas_pos[i] - i
            rho = self._collapse_single(rho, pos,
                                        projections[i], "Z")
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

    def _measure_random_ancillas_Z(self, rho, ancillas_pos):
        # NOTE: All ancillas are collapsed in parallel
        self.check["measurement"] += 1

        time = self.time_lookup["measurement"]
        self.check["time"] += time
        self.check["time1"] += time
        # Apply environmental error
        rho = errs.env_error_all(rho, 0, self.a1, time)

        # Collapse the qubits in parrallel
        # Sort list to be able to reduce dimension and keep track of positions
        ancillas_pos = sorted(ancillas_pos)
        measurements = []
        for i in range(len(ancillas_pos)):
            pos = ancillas_pos[i] - i
            m, rho = self._measure_random_single(rho, pos, "Z")
            measurements += [m]
        return measurements, rho

    def _measure_random_single(self, rho, pos, basis):
        # Measure a single qubit in the state.
        # NOTE: Dimension is reduced after collapse
        N = len(rho.dims[0])

        if basis == "X":
            rho = errs.single_qubit_gate_noise(rho, self.ps, N, pos)
            measurement, rho = errs.measure_single_Xbasis_random(rho, self.pm,
                                                                 N, pos)
        elif basis == "Z":
            measurement, rho = errs.measure_single_Zbasis_random(rho, self.pm,
                                                                 N, pos)
        return measurement, rho


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

    def _apply_single_qubit_gates(self, rho, gates, operation_qubits):
        # Apply single qubit gates on the entire state rho
        N = len(rho.dims[0])
        # NOTE: Number of gates is 1 because gates are applied in parallel
        # in the nodes
        self.check["single_qubit_gate"] += 1
        time = self.time_lookup["single_qubit_gate"]
        self.check["time"] += time
        self.check["time1"] += time
        # Apply environmental error
        rho = errs.env_error_all(rho, 0, self.a1, time)

        for i in range(len(gates)):
            rho = errs.single_qubit_gate(rho, gates[i], self.ps, N,
                                         operation_qubits[i])
        return rho

    def start_bell_pair(self, rho=None):
        """
        Start with a raw Bell pair using the single click protocol.

        NOTE: rho is not taken from arguments, only exists due to
        how 'circuit.py' is constructed.
        """
        self._reset_check()
        time, bell = self._generate_bell_single_click()
        self.check["time"] += time
        self.check["time0"] += time
        return 1, self.check, bell

    def start_BK(self, rho=None):
        """
        Circuit that generates a Bell pair from the Barret-Kok protocol.

        NOTE: rho is not taken from arguments, only exists due to
        how 'circuit.py' is constructed.
        """
        self._reset_check()
        time, bell = self._generate_bell_pair_BK()
        self.check["time"] += time
        self.check["time0"] += time
        return 1, self.check, bell

    def start_epl(self, rho=None):
        """
        Circuit that generates a Bell pair using the EPL protocol.

        NOTE: rho is not taken from arguments, only exists due to
        how 'circuit.py' is constructed.
        """
        self._reset_check()
        # Generate Bell pairs
        # First pair with swap
        time1, bell1 = self._generate_bell_single_click()
        bell1 = self._swap_pair(bell1, [0, 1])
        # Second pair
        time2, bell2 = self._generate_bell_single_click()
        self.check["time"] += time1 + time2
        self.check["time0"] += time1 + time2
        # Apply cutoff here

        # Dephase pair 1
        rho = errs.env_error_all(bell1, self.a0, self.a1, time2)

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
        p_success, rho = self._collapse_ancillas_EPL(rho, targets)

        # Apply X to rotate state
        rho = self._apply_single_qubit_gates(rho, [qt.rx(np.pi, 2, 0)], [0])

        return p_success, self.check, rho

    def add_bell_pair(self, rho):
        """Append a single click Bell pair to the state.

        Parameters
        ----------
        rho : (densmat) density matrix in which the bell pair appends.
        """
        self._reset_check()
        time, bell = self._generate_bell_single_click()
        # Apply environmental error
        rho = errs.env_error_all(rho, self.a0, self.a1, time)
        self.check["time"] += time
        self.check["time1"] += time
        # Bell state is attached at the end of the complete state
        rho = qt.tensor(rho, bell)
        return 1, self.check, rho

    def swap_pair(self, rho, pair):
        """
        Does not really make a SWAP, only applies the noise.
        Position of states in each node must be keep track manually.

        Parameters
        ----------
        rho : (densmat) density matrix
        pair : (list) pair of qubits to be swaped ex. [2, 3]
        """
        self._reset_check()
        # Apply the noise
        self._swap_pair(rho, pair)
        return 1, self.check, rho

    def two_qubit_gates(self, rho, controls, targets, sigma):
        """
        Apply one local two qubit Control type of gates in parallel
        on each node.

        Parameters
        ----------
        rho : (densmat) density matrix
        controls : (list) list with the control qubits position
        targets : (list) list with the target qubits position
        sigma : (string) X or Z depending on the gate
        """
        self._reset_check()
        rho = self._apply_two_qubit_gates(rho, controls, targets, sigma)
        return 1, self.check, rho

    def collapse_ancillas_X(self, rho, ancillas_pos, projections):
        """
        Measure the ancillas in the X basis in parallel in each node.
        Ancillas position need to be the last part of the state.

        Parameters
        ----------
        rho : (densmat) density matrix
        ancillas_pos : (int) list with the positions of ancillas to be measured
        projections : (list) list with the projections in which every ancilla collapses
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

        Parameters
        ----------
        rho : (densmat) density matrix
        ancillas_pos : (int) list with the positions of ancillas to be measured
        projections : (list) list with the projections in which every ancilla collapses
        """
        # Reset check
        self._reset_check()
        p_success, rho = self._collapse_ancillas_Z(rho, ancillas_pos,
                                                   projections)
        return p_success, self.check, rho

    def collapse_ancillas_GHZ(self, rho, ghz_size, measure_pos):
        """
        Collapse the ancillas in the nodes to create a GHZ state.
        In the ghz_4 connected through 4 pairs then the
        ancillas must be the last qubit(s)

        Parameters
        ----------
        rho : (densmat) density matrix
        ghz_size : (int) size of the GHZ state being formed, must be 3 or 4
        measure_pos : (int) list with the positions of ancillas to be measured
        """
        # Reset number of steps counter
        self._reset_check()

        N = len(rho.dims[0])
        p_success = 1

        # Special case: Success probability for the GHZ4 case with 4 pairs.
        if ghz_size == 4 and len(measure_pos) == 4:
            p_success = ps.ghz_4(rho)

        # Get random measuements outcomes
        measurements, rho_measured = self._measure_random_ancillas_Z(rho, measure_pos)
        # Transform measurements from 1 and -1 to 0 and 1
        measurements = np.array(measurements) % 3 - 1
        # print(measurements)
        N = len(rho_measured.dims[0])
        # The qubits in which the correction applies
        if ghz_size == 3:
            # The last qubit if GHZ3
            operation_pos = N - 1
        elif ghz_size == 4:
            # The last two if GHZ4
            operation_pos = [N - 1, N - 2]

        correction = gc.correction(measurements, ghz_size, N, operation_pos)
        # For the case of GHZ4 using 4 pairs theres a chance the entire protocol
        # fails due the to the extra error correction, if thats the case just
        # redo measurements.
        while correction is None:
            measurements, rho_measured = self._measure_random_ancillas_Z(rho, measure_pos)
            measurements = np.array(measurements) % 3 - 1
            correction = gc.correction(measurements, ghz_size, N, operation_pos)

        rho = rho_measured
        # Check if correction is really nedded
        if len(correction) != 0:
            rho = self._apply_single_qubit_gates(rho, correction, operation_pos)

        # print(p_success)
        # print(measurements)
        return p_success, self.check, rho

    def single_selection(self, rho, operation_qubits, sigma):
        """
        Single selection round. Uses 2 ancilla qubits.

        Parameters
        ----------
        rho : (densmat) density matrix
        operation_pos : (list) list with the positions of the qubits to be purified
                        by the protocol
        sigma : (string) X or Z parity used in the purification
        """
        # Reset number of steps counter
        self._reset_check()

        # Generate raw bell pair
        rho = self._append_epl(rho)
        N = len(rho.dims[0])

        H = [qt.snot(N, operation_qubits[0]), qt.snot(N, operation_qubits[1])]
        rho = self._apply_single_qubit_gates(rho, H, operation_qubits)

        # Apply two qubit gates
        controls = [N-1, N-2]
        # self._swap_pair(rho, controls)
        # rho = self._apply_two_qubit_gates(rho, controls,
        #                                   operation_qubits, sigma)
        rho = self._apply_two_qubit_gates(rho, operation_qubits,
                                          controls, sigma)

        # Calculate the probability of success
        p_success = ps.single_sel(rho, N, controls)
        # Measure ancillas in X basis
        projections = [0, 0]
        _, rho = self._collapse_ancillas_X(rho, controls, projections)
        return p_success, self.check, rho

    def double_selection(self, rho, operation_qubits, sigma):
        """
        Double selection round. Uses 4 ancilla qubits.

        Parameters
        ----------
        rho : (densmat) density matrix
        operation_pos : (list) list with the positions of the qubits to be purified
                        by the protocol
        sigma : (string) X or Z parity used in the purification
        """
        # Reset number of steps counter
        self._reset_check()

        # Generate two bell pairs with a SWAP in between them
        rho = self._append_epl(rho)
        N = len(rho.dims[0])
        self._swap_pair(rho, [N-1, N-2])
        rho = self._append_epl(rho)
        N += 2

        # Apply first two qubit gates
        controls1 = [N-3, N-4]
        rho = self._apply_two_qubit_gates(rho, controls1,
                                          operation_qubits, sigma)

        # Apply second set of gates
        controls2 = [N-1, N-2]
        # Swap noise
        self._swap_pair(rho, controls2)
        rho = self._apply_two_qubit_gates(rho, controls2, controls1, "Z")

        # Success probability
        p_success = ps.double_sel(rho, N, controls1, controls2)

        # Measure ancillas in X basis
        projections = [0] * 2
        p_success2, rho = self._collapse_ancillas_X(rho, controls2, [0, 0])
        # Swap noise
        self._swap_pair(rho, controls1)
        p_success1, rho = self._collapse_ancillas_X(rho, controls1, projections)

        return p_success, self.check, rho


    def double_selection22(self, sigma):
        """
        Double selection round. Uses 4 ancilla qubits.

        Parameters
        ----------
        rho : (densmat) density matrix
        operation_pos : (list) list with the positions of the qubits to be purified
                        by the protocol
        sigma : (string) X or Z parity used in the purification
        """
        # Reset number of steps counter
        self._reset_check()

        # Generate raw bell pair
        _, _, rho = self.start_epl()
        rho = self._append_epl(rho)
        N = len(rho.dims[0])
        operation_qubits = [0, 1]

        H = [qt.snot(N, operation_qubits[0]), qt.snot(N, operation_qubits[1])]
        rho = self._apply_single_qubit_gates(rho, H, operation_qubits)

        # Apply two qubit gates
        controls = [N-1, N-2]
        # self._swap_pair(rho, controls)
        # rho = self._apply_two_qubit_gates(rho, controls,
        #                                   operation_qubits, sigma)
        rho = self._apply_two_qubit_gates(rho, operation_qubits,
                                          controls, "Z")

        # Calculate the probability of success
        p_success1 = ps.single_sel(rho, N, controls)
        # Measure ancillas in X basis
        projections = [0, 0]
        _, rho = self._collapse_ancillas_X(rho, controls, projections)

        rho = self._append_epl(rho)
        N = len(rho.dims[0])

        H = [qt.snot(N, operation_qubits[0]), qt.snot(N, operation_qubits[1])]
        rho = self._apply_single_qubit_gates(rho, H, operation_qubits)

        # Apply two qubit gates
        controls = [N-1, N-2]
        # self._swap_pair(rho, controls)
        # rho = self._apply_two_qubit_gates(rho, controls,
        #                                   operation_qubits, sigma)
        rho = self._apply_two_qubit_gates(rho, operation_qubits,
                                          controls, sigma)

        # Calculate the probability of success
        p_success2 = ps.single_sel(rho, N, controls)
        # Measure ancillas in X basis
        projections = [0, 0]
        _, rho = self._collapse_ancillas_X(rho, controls, projections)

        p_success = p_success1 * p_success2

        return p_success, self.check, rho
