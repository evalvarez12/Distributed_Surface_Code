"""
Noise modeling using the Choi-Jamiolkowski isomorphism
to find the chi matrix that caracterizes the process

created-on: 30/06/17
"""
import numpy as np
import qutip as qt
import itertools
import pauli_basis
import operations as ops

class NoiseModel:
    """
    Class to hold all the thigs used.
    """

    def __init__(self, system_size, superoperator_function, parity):
        """Init fucntion to prepare all required assets."""
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
        self.rho = self.psi_basis[0] * self.psi_basis[0].dag()
        self.rho = superoperator_function(self.rho, self.targets, self.basis_parity)

        # Empty dictionary to store the chi matrix
        self.chi = {}

    def _remove_sym_pauli_basis(self):
        sym = self.basis_parity * self.system_size

        syms = list(self.pauli_basis.keys())
        for k in syms:
            if k not in self.pauli_basis:
                continue

            k_equiv = pauli_basis.symbol_product(k, sym)

            if k_equiv not in self.pauli_basis:
                continue

            # To improve readability keep the simpliest symbol
            if "I" in k:
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
                symOK = k + "OK"
                self.chi[symOK] = self._fidelity(self.rho, v * self.psi_basis[0])

                # Decomposition when there is a faulty measurement
                symNOK = k + "NOK"
                self.chi[symNOK] = self._fidelity(self.rho, v * self.psi_basis[1])
            else:
                # Decompose through all the Pauli basis
                self.chi[k] = self._fidelity(self.rho, v * self.psi_basis[0])
        # Call function to reduce permuations
        self._chi_reduce_permutations()

    def check_total_num(self):
        total = 0
        for v in self.chi.values():
            total += v
        return total

    def _chi_reduce_permutations(self):
        syms = list(self.chi.keys())
        for k in syms:
            if k not in self.chi:
                continue

            # print("KEY: ", k)
            symbol = k
            p_sum = 0
            residue = ""

            if self.faulty_measurement:
                if "N" in symbol:
                    symbol = symbol.replace("NOK", "")
                    residue = "NOK"
                else:
                    symbol = symbol.replace("OK", "")
                    residue = "OK"

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
                # If success try with the parity correspondant
                # else:
                #     str_symbol = pauli_basis.symbol_product(str_symbol,
                #                                             self.basis_parity)
                #     p_symbol = str_symbol + residue
                #     print("Will try thi one: ", p_symbol)
                #     p_sum += self.chi[p_symbol]
                #     del self.chi[p_symbol]

            if p_sum != 0:
                self.chi[k] = p_sum














################################################# CODE LIMBO #############################################

#         for err_type in pauli_errs:
#             # Save the names of the errors
#             type_name = err_type[1]
#
#             # Array to store all fidelities results
#             weights = []
#
#             # Odd and even reference collapsed states
#             even_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)
#             odd_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)
#
#             for permut in err_type[0]:
#                 # errs_permut = []
#                 err_permut_good = 0
#                 err_permut_lie = 0
#
#                 for e in permut:
#                     # e_psi0 = e[0] * even_psi
#                     # e_psi1 = e[0] * odd_psi
#                     e_psi0 = e[0] * psi_initial
#                     e_psi1 = e[0] * psi_initial
#
#
#                     # print(prot.fidelity(rhos[0], e_psi0))
#                     err_permut_good += (prot.fidelity(rhos[0], e_psi0))
#                                         # + ps[1]*prot.fidelity(rhos[1], e_psi1)
#
#                     # err_permut_lie += (ps[0]*prot.fidelity(rhos[0], e_psi1)
#                     #                    + ps[1]*prot.fidelity(rhos[1], e_psi0))
#
#                     err_permut_lie += (0)
#
#                 weights += [err_permut_good, err_permut_lie, e[1]]
#
#             # Convert to arrays
#             weights_array = np.array(weights)
#
#             # Save arrays
#             np.save("data/decomposition_" + type_name, weights_array)
#
#
# #
# # # Set system size
# # system_size = 4
# #
# # # Initialize protocols with errors
# # prot = protocols.Protocols(ps=0.0, pm=0.009, pg=0.009, pn=0.)
# #
# # # The initial state
# # psi_initial = prot.operational_state_ket(system_size)
# # rho_initial = psi_initial * psi_initial.dag()
# #
# # # Parity measurement is made on one side of the bipartite state
# # targets = [0, 2, 4, 6]
# #
# # parity = "Z"
# # # Get the probabilities and the collapsed states
# # ps_seq, rhos = prot.monolithic(rho_initial, targets, parity)
# # ps = np.average(ps_seq, axis=1)
# #
# #
# # # Get all errors
# # pauli_errs = pauli_errors.get_pauli_errors(targets, 2*system_size)
# #
# # # Funtion to realize the decomposition
# # def decompose(ps, rhos, psi_initial):
# #     for err_type in pauli_errs:
# #         # Save the names of the errors
# #         type_name = err_type[1]
# #
# #         # Array to store all fidelities results
# #         weights = []
# #
# #         # Odd and even reference collapsed states
# #         even_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)
# #         odd_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)
# #
# #         for permut in err_type[0]:
# #             # errs_permut = []
# #             err_permut_good = 0
#             err_permut_lie = 0
#
#             for e in permut:
#                 # e_psi0 = e[0] * even_psi
#                 # e_psi1 = e[0] * odd_psi
#                 e_psi0 = e[0] * psi_initial
#                 e_psi1 = e[0] * psi_initial
#
#
#                 # print(prot.fidelity(rhos[0], e_psi0))
#                 err_permut_good += (prot.fidelity(rhos[0], e_psi0))
#                                     # + ps[1]*prot.fidelity(rhos[1], e_psi1)
#
#                 # err_permut_lie += (ps[0]*prot.fidelity(rhos[0], e_psi1)
#                 #                    + ps[1]*prot.fidelity(rhos[1], e_psi0))
#
#                 err_permut_lie += (0)
#
#             weights += [err_permut_good, err_permut_lie, e[1]]
#
#         # Convert to arrays
#         weights_array = np.array(weights)
#
#         # Save arrays
#         np.save("data/decomposition_" + type_name, weights_array)
#
#
# decompose(ps, rhos, psi_initial)
# # decompose(ps, [rho_initial]*2, psi_initial)
#
#
# # Load all data
# decomposeI = np.load("data/decomposition_I.npy")
# decomposeA = np.load("data/decomposition_A.npy")
# decomposeB = np.load("data/decomposition_B.npy")
# decomposeC = np.load("data/decomposition_C.npy")
# decomposeD = np.load("data/decomposition_D.npy")
#
#
# # Remove the string names that occur every third element
# PdecomposeI = decomposeI[np.mod(np.arange(len(decomposeI)), 3) != 2]
# PdecomposeA = decomposeA[np.mod(np.arange(len(decomposeA)), 3) != 2]
# PdecomposeB = decomposeB[np.mod(np.arange(len(decomposeB)), 3) != 2]
# PdecomposeC = decomposeC[np.mod(np.arange(len(decomposeC)), 3) != 2]
# PdecomposeD = decomposeD[np.mod(np.arange(len(decomposeD)), 3) != 2]
#
#
# # The ones with D have a systematic err that needs to be removed,(XXXX or ZZZZ)
# # PdecomposeD = PdecomposeD[PdecomposeD != PdecomposeI[0]]
# # PdecomposeD = PdecomposeD[PdecomposeD != PdecomposeI[1]]
#
# GOODdecomposeI = PdecomposeI[::2].astype(np.float64)
# GOODdecomposeA = PdecomposeA[::2].astype(np.float64)
# GOODdecomposeB = PdecomposeB[::2].astype(np.float64)
# GOODdecomposeC = PdecomposeC[::2].astype(np.float64)
# GOODdecomposeD = PdecomposeD[::2].astype(np.float64)
#
#
# # Convert arrays to floats
# f_decomposeI = PdecomposeI.astype(np.float64)
# f_decomposeA = PdecomposeA.astype(np.float64)
# f_decomposeB = PdecomposeB.astype(np.float64)
# f_decomposeC = PdecomposeC.astype(np.float64)
# f_decomposeD = PdecomposeD.astype(np.float64)
#
# total = (f_decomposeA.sum() + f_decomposeB.sum() + f_decomposeC.sum()
#          + f_decomposeD.sum() + f_decomposeI.sum())
#
#
# totalGOOD = (GOODdecomposeA.sum() + GOODdecomposeB.sum() + GOODdecomposeC.sum()
#              + GOODdecomposeD.sum() + GOODdecomposeI.sum())
#
# # badD = badD[badD != badD.max()]
# #
# # goods = 0
# # for i in [goodA, goodB, goodC, goodD]:
# #     goods += i.sum()
# #
# # bads = 0
# # for i in [badA, badB, badC, badD]:
# #     bads += i.sum()
