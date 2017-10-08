"""
Noise modeling using the Choi-Jamiolkowski isomorphism
to find the chi matrix that caracterizes the process

created-on: 30/06/17
"""
import numpy as np
import qutip as qt
import pickle
import itertools
from . import pauli_basis
from . import operations as ops

class NoiseModel:
    """
    Class to hold all the thigs used.
    """

    def __init__(self, system_size, superoperator_function, parity):
        """Init function to prepare all required assets."""
        # Create the initial state
        self.psi_basis = [self._choi_state_ket(system_size)]
        self.faulty_measurement = False
        self.basis_parity = parity
        # The state comes organized in an ordered manner
        self.targets = list(range(system_size))
        self.pauli_basis = pauli_basis.get_basis(self.targets, 2 * system_size)

        # Save initial variables
        self.system_size = system_size

        # Apply the superoperator to the choi state
        self.rhos = self.psi_basis[0] * self.psi_basis[0].dag()
        self.ps, self.rhos = superoperator_function(self.rhos, self.targets, self.basis_parity)

        # Empty dictionary to store the chi matrix
        self.chi = {}

    def save_chi(self, file_name):
        pickle_out = open("data/" + file_name, "wb")
        pickle.dump(self.chi, pickle_out)
        pickle_out.close()

    def _remove_sym_pauli_basis(self):
        sym = self.basis_parity * self.system_size

        syms = list(self.pauli_basis.keys())
        for k in syms:
            if k not in self.pauli_basis:
                continue

            k_equiv = pauli_basis.symbol_product(k, sym)

            if k_equiv not in self.pauli_basis:
                continue

            # To improve readability keep the simpNOKst symbol
            if k.count("I") > k_equiv.count("I"):
                del self.pauli_basis[k_equiv]
            else:
                del self.pauli_basis[k]

    def _fidelity(self, rhoA, stateB):
        """
        Fidelity for the special case when one of the states is a pure state.
        """
        return (stateB.dag() * rhoA * stateB).norm()

    def _parity_projection_ket(self, psi, targets, measurement):
        plus = qt.snot() * qt.basis(2, 0)
        full_psi = qt.tensor(psi, plus)
        N = len(full_psi.dims[0])
        control = N-1
        if self.basis_parity == "X":
            for t in targets:
                full_psi = qt.cnot(N, control, t) * full_psi
        if self.basis_parity == "Z":
            for t in targets:
                full_psi = qt.cphase(np.pi, N, control, t) * full_psi

        # NOTE: Parity projection is not a channel
        collapsed_psi = ops.collapse_single_Xbasis_ket(full_psi, measurement,
                                                       N, N-1, True)
        if collapsed_psi.norm() != 0:
            collapsed_psi = collapsed_psi/collapsed_psi.norm()
        return collapsed_psi

    def separate_basis_parity(self):
        psi_even = self._parity_projection_ket(self.psi_basis[0],
                                               self.targets,
                                               0)
        psi_odd = self._parity_projection_ket(self.psi_basis[0],
                                              self.targets,
                                              1)

        # Separate into even and odd parity
        self.psi_basis = [psi_even, psi_odd]
        # Remove the duplicated pauli vecs
        self._remove_sym_pauli_basis()
        # Mark that parity symmetry is used to indicate faulty_measurement
        self.faulty_measurement = True

    def _choi_state_ket(self, N):
        """
        Create a bipartite state.
        Used in the C-J isomorphism.
        """
        state = qt.bell_state('00')
        for i in range(1, N):
            state = qt.tensor(qt.bell_state('00'), state)
        # The desired order of the qubits
        order = list(range(0, 2*N, 2)) + list(range(1, 2*N, 2))
        state = state.permute(order)
        return state

    def make_chi_matrix(self):
        for k, v in self.pauli_basis.items():
            if self.faulty_measurement:
                # Decomposition when no faulty measurement
                sym_ok = k + "_OK"
                f_even_ok = self._fidelity(self.rhos[0], v * self.psi_basis[0])
                f_odd_ok = self._fidelity(self.rhos[1], v * self.psi_basis[1])
                self.chi[sym_ok] = (self.ps[0] * f_even_ok
                                   + self.ps[1] * f_odd_ok)
                # Decomposition when there is a faulty measurement
                sym_NOK = k + "_NOK"
                f_even_NOK = self._fidelity(self.rhos[0], v * self.psi_basis[1])
                f_odd_NOK = self._fidelity(self.rhos[1], v * self.psi_basis[0])
                self.chi[sym_NOK] = (self.ps[0] * f_even_NOK
                                    + self.ps[1] * f_odd_NOK)
            else:
                # Decompose through all the Pauli basis
                self.chi[k] = self._fidelity(self.rho, v * self.psi_basis[0])
        # Call function to reduce permuations
        self._chi_reduce_permutations()

    def check_total_sum(self):
        total = 0
        for v in self.chi.values():
            total += v
        return total

    def reset_chi(self):
        self.chi = {}

    def _chi_reduce_permutations(self):
        syms = list(self.chi.keys())
        prod_parity = self.basis_parity * self.system_size
        for k in syms:
            if k not in self.chi:
                continue

            # print("KEY: ", k)
            symbol = k
            p_sum = 0
            residue = ""

            if self.faulty_measurement:
                if "N" in symbol:
                    symbol = symbol.replace("_NOK", "")
                    residue = "_NOK"
                else:
                    symbol = symbol.replace("_OK", "")
                    residue = "_OK"

            perms = [''.join(p) for p in itertools.permutations(symbol)]
            perms = set(perms)
            for p in perms:
                str_symbol = ''.join(p)

                p_symbol = str_symbol + residue
                # print("Peruting: ", p_symbol)
                # Try with the current symbol permutation
                if p_symbol in self.chi:
                    # print("He is in chi: ", p_symbol)
                    p_sum += self.chi[p_symbol]
                    del self.chi[p_symbol]
                # If not success try with the parity correspondant
                else:
                    str_symbol = pauli_basis.symbol_product(str_symbol,
                                                            prod_parity)
                    p_symbol = str_symbol + residue
                    # print("Will try thi one: ", p_symbol)
                    if p_symbol in self.chi:
                        p_sum += self.chi[p_symbol]
                        del self.chi[p_symbol]

            if p_sum != 0:
                self.chi[k] = p_sum
