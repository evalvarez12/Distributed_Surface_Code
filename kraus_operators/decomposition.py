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
ps_seq, rhos = prot.monolithic(rho_initial, targets, parity)
ps = np.average(ps_seq, axis=1)

# Odd and even reference collapsed states
even_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)
odd_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)

# Get all errors
pauli_errs = pauli_errors.get_pauli_errors(targets, 2*system_size)

# Funtion to realize the decomposition
def decompose(ps, rhos, psis):
    for err_type in pauli_errs:
        # Save the names of the errors
        type_name = err_type[1]

        # Array to store all fidelities results
        weights = []

        for permut in err_type[0]:
            # errs_permut = []
            err_permut_good = 0
            err_permut_lie = 0

            for e in permut:
                e_psi0 = e[0] * psis[0]
                e_psi1 = e[0] * psis[1]

                # print(prot.fidelity(rhos[0], e_psi0))
                err_permut_good += (ps[0]*prot.fidelity(rhos[0], e_psi0)
                                    + ps[1]*prot.fidelity(rhos[1], e_psi1))

                err_permut_lie += (ps[0]*prot.fidelity(rhos[0], e_psi1)
                                   + ps[1]*prot.fidelity(rhos[1], e_psi0))

            weights += [err_permut_good, err_permut_lie, e[1]]

        # Convert to arrays
        weights_array = np.array(weights)

        # Save arrays
        np.save("data/decomposition_" + type_name, weights_array)


decompose(ps, rhos, [even_psi, odd_psi])

# Load all data
decomposeI = np.load("data/decomposition_I.npy")
decomposeA = np.load("data/decomposition_A.npy")
decomposeB = np.load("data/decomposition_B.npy")
decomposeC = np.load("data/decomposition_C.npy")
decomposeD = np.load("data/decomposition_D.npy")


# Remove the string names that occur every third element
PdecomposeI = decomposeI[np.mod(np.arange(len(decomposeI)), 3) != 2]
PdecomposeA = decomposeA[np.mod(np.arange(len(decomposeA)), 3) != 2]
PdecomposeB = decomposeB[np.mod(np.arange(len(decomposeB)), 3) != 2]
PdecomposeC = decomposeC[np.mod(np.arange(len(decomposeC)), 3) != 2]
PdecomposeD = decomposeD[np.mod(np.arange(len(decomposeD)), 3) != 2]


# The ones with D have a systematic err that needs to be removed,(XXXX or ZZZZ)
# PdecomposeD = PdecomposeD[PdecomposeD != PdecomposeI[0]]
# PdecomposeD = PdecomposeD[PdecomposeD != PdecomposeI[1]]

GOODdecomposeI = PdecomposeI[::2].astype(np.float64)
GOODdecomposeA = PdecomposeA[::2].astype(np.float64)
GOODdecomposeB = PdecomposeB[::2].astype(np.float64)
GOODdecomposeC = PdecomposeC[::2].astype(np.float64)
GOODdecomposeD = PdecomposeD[::2].astype(np.float64)


# Convert arrays to floats
f_decomposeI = PdecomposeI.astype(np.float64)
f_decomposeA = PdecomposeA.astype(np.float64)
f_decomposeB = PdecomposeB.astype(np.float64)
f_decomposeC = PdecomposeC.astype(np.float64)
f_decomposeD = PdecomposeD.astype(np.float64)

total = (f_decomposeA.sum() + f_decomposeB.sum() + f_decomposeC.sum()
         + f_decomposeD.sum() + f_decomposeI.sum())


totalGOOD = (GOODdecomposeA.sum() + GOODdecomposeB.sum() + GOODdecomposeC.sum()
             + GOODdecomposeD.sum() + GOODdecomposeI.sum())

# badD = badD[badD != badD.max()]
#
# goods = 0
# for i in [goodA, goodB, goodC, goodD]:
#     goods += i.sum()
#
# bads = 0
# for i in [badA, badB, badC, badD]:
#     bads += i.sum()
