"""
Noise modeling using the Choi-Jamiolkowski isomorphism
to find the chi matrix that caracterizes the process

created-on: 30/06/17
"""
import numpy as np
import qutip as qt
import protocols
import pauli_basis
import projectors


class NoiseModel:
    """
    Class to hold all the thigs used.
    """

    def __init__(self, system_size, superoperator_function, symmetry):
        """Init fucntion to prepare all required assets."""
        # Create the initial state
        self.psi_initial = self.choi_state_ket(system_size)
        # The state comes organized in an ordered manner
        self.targets = list(range(system_size))
        self.rho = self.psi_initial * self.psi_initial.dag()
        self.pauli_basis = pauli_basis.get_basis(self.targets, system_size)

        # Save initial variables
        self.system_size = system_size
        self.superoperator = superoperator_function

        # Empty dictionary to store the chi matrix
        self.model = {}

    def remove_sym_pauli_basis(self, symmetry):
        if symmetry = "X":
            sym = "X" * self.system_size
        if symmetry = "Z":
            sym = "Z" * self.system_size

        for k in self.pauli_basis.keys():
            k_prod = pauli_basis.symbol_product(k, sym)

            # To improve readability keep the simpliest symbol
            if "I" in k:
                del self.pauli_basis[k_prod]
            else:
                del self.pauli_basis[k]

    def fidelity(self, rhoA, stateB):
        """
        Fidelity for the special case when one of the states is a pure state.
        """
        return (stateB.dag() * rhoA * stateB).norm()

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

        # NOTE: Parity projection is not a channel
        collapsed_psi = ops.collapse_single_Xbasis_ket(full_psi, measurement,
                                                       N, N-1, True)
        return collapsed_psi

    def choi_state_ket(self, N):
        """
        Create a bipartite state.
        Used in the C-J isomorphism.
        """
        state = qt.bell_state('00')
        for i in range(1, N):
            state = qt.tensor(qt.bell_state('00'), state)
        # The desired order of the qubits
        order = list(range(0, 2*N, 2)) + list(range(1, 2*N, 2))
        state = state.permut(order)
        return state

    def make_chi_matrix(self):
        for err_type in pauli_errs:
            # Save the names of the errors
            type_name = err_type[1]

            # Array to store all fidelities results
            weights = []

            # Odd and even reference collapsed states
            even_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)
            odd_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)

            for permut in err_type[0]:
                # errs_permut = []
                err_permut_good = 0
                err_permut_lie = 0

                for e in permut:
                    # e_psi0 = e[0] * even_psi
                    # e_psi1 = e[0] * odd_psi
                    e_psi0 = e[0] * psi_initial
                    e_psi1 = e[0] * psi_initial


                    # print(prot.fidelity(rhos[0], e_psi0))
                    err_permut_good += (prot.fidelity(rhos[0], e_psi0))
                                        # + ps[1]*prot.fidelity(rhos[1], e_psi1)

                    # err_permut_lie += (ps[0]*prot.fidelity(rhos[0], e_psi1)
                    #                    + ps[1]*prot.fidelity(rhos[1], e_psi0))

                    err_permut_lie += (0)

                weights += [err_permut_good, err_permut_lie, e[1]]

            # Convert to arrays
            weights_array = np.array(weights)

            # Save arrays
            np.save("data/decomposition_" + type_name, weights_array)


#
# # Set system size
# system_size = 4
#
# # Initialize protocols with errors
# prot = protocols.Protocols(ps=0.0, pm=0.009, pg=0.009, pn=0.)
#
# # The initial state
# psi_initial = prot.operational_state_ket(system_size)
# rho_initial = psi_initial * psi_initial.dag()
#
# # Parity measurement is made on one side of the bipartite state
# targets = [0, 2, 4, 6]
#
# parity = "Z"
# # Get the probabilities and the collapsed states
# ps_seq, rhos = prot.monolithic(rho_initial, targets, parity)
# ps = np.average(ps_seq, axis=1)
#
#
# # Get all errors
# pauli_errs = pauli_errors.get_pauli_errors(targets, 2*system_size)
#
# # Funtion to realize the decomposition
# def decompose(ps, rhos, psi_initial):
#     for err_type in pauli_errs:
#         # Save the names of the errors
#         type_name = err_type[1]
#
#         # Array to store all fidelities results
#         weights = []
#
#         # Odd and even reference collapsed states
#         even_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)
#         odd_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)
#
#         for permut in err_type[0]:
#             # errs_permut = []
#             err_permut_good = 0
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
