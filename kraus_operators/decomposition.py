"""
Decomposition routines into Kraus operators.

created-on: 30/06/17
"""
import numpy as np
import protocols
import pauli_errors
import projectors


# Set system size
system_size = 4

# Initialize protocols with errors
prot = protocols.Protocols(ps=0.0, pm=0.009, pg=0.009, pn=0.)

# The initial state
psi_initial = prot.operational_state_ket(system_size)
rho_initial = psi_initial * psi_initial.dag()

# Parity measurement is made on one side of the bipartite state
targets = [0, 2, 4, 6]

parity = "Z"
# Get the probabilities and the collapsed states
ps, rhos = prot.monolithic(rho_initial, targets, parity)
# print("measurement: ", measurement)

# Odd and even reference collapsed states
even_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)
odd_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)



f1 = prot.fidelity(rhos[0], even_psi)
f2 = prot.fidelity(rhos[1], odd_psi)
print(f1)
print(f2)
print(ps[0][0]*f1 + ps[1][0]*f2)
print("LIE")
print(prot.fidelity(rhos[1], even_psi))
print(prot.fidelity(rhos[0], odd_psi))

# Get all errors
pauli_errs = pauli_errors.get_pauli_errors(2*system_size, targets)

# Funtion to realize the decomposition
# def decompose(rhos, psis, subfix):
#     for err_type in pauli_errs:
#         # Save the names of the errors
#         names = err_type[1][1]
#
#         # Array to store all fidelities results
#         errs = []
#
#         # Do the comparation with no pauli error
#         if err_type[0] == "A":
#             errs += [prot.fidelity(rhos[0], psis[0])] + [prot.fidelity(rhos[1], psis[1])]
#             names = ["I"] + names
#
#         for permut in err_type[1][0]:
#             # errs_permut = []
#             errs_permut = 0
#             for e in permut:
#                 e_psi0 = e * psis[0]
#                 e_psi1 = e * psis[1]
#
#                 # errs_permut += [prot.fidelity(rho, e_psi)]
#                 errs_permut += prot.fidelity(rhos[0], e_psi0) + prot.fidelity(rhos[1], e_psi1)
#             errs += [errs_permut]
#
#         # Convert to arrays
#         errs_array = np.array(errs)
#         names_array = np.array(names)
#
#         # Save arrays
#         np.save("data/errors_" + subfix + err_type[0], errs_array)
#         np.save("data/errors_names_" + subfix + err_type[0], names_array)

#
# decompose(rhos, [good_psi, bad_psi], "good_")
# # decompose(rhos, bad_psi, measurement, "bad_")
#
#
# # Load all data
# goodA = np.load("data/errors_good_A.npy")
# badA = np.load("data/errors_bad_A.npy")
# goodB = np.load("data/errors_good_B.npy")
# badB = np.load("data/errors_bad_B.npy")
# goodC = np.load("data/errors_good_C.npy")
# badC = np.load("data/errors_bad_C.npy")
# goodD = np.load("data/errors_good_D.npy")
# badD = np.load("data/errors_bad_D.npy")
#
#
# # The ones with D have a systematic err that needs to be removed
# goodD = goodD[goodD != goodA.max()]
# badD = badD[badD != badD.max()]
#
# goods = 0
# for i in [goodA, goodB, goodC, goodD]:
#     goods += i.sum()
#
# bads = 0
# for i in [badA, badB, badC, badD]:
#     bads += i.sum()
