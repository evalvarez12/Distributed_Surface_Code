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
prot = protocols.Protocols(ps=0.0, pm=0.006, pg=0.006, pn=0.1)

# The initial state
psi_initial = prot.operational_state_ket(system_size)
rho_initial = psi_initial * psi_initial.dag()

# Parity measurement is made on one side of the bipartite state
targets = [0, 2, 4, 6]

parity = "X"
measurement, rho = prot.expedient(rho_initial, targets, parity)
print("measurement: ", measurement)




# Measurements to array to use .prod()
measurement = np.array(measurement).prod()

# Do the parity distintion to identify the good vs bad measurement
if measurement.prod() == 1:
    good_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)
    bad_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)
if measurement.prod() == -1:
    good_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)
    bad_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)


# Get all errors
pauli_errs = pauli_errors.get_pauli_errors(2*system_size, targets)

# Funtion to realize the decomposition
def decompose(rho, psi, measurement, subfix):
    for err_type in pauli_errs:
        # Save the names of the errors
        names = err_type[1][1]

        # Array to store all fidelities results
        errs = []

        # Do the comparation with no pauli error
        if err_type[0] == "A":
            errs += [prot.fidelity(rho, psi)]
            names = ["I"] + names

        for permut in err_type[1][0]:
            # errs_permut = []
            errs_permut = 0
            for e in permut:
                e_psi = e * psi

                # # Do the parity distintion to identify the good vs bad measurement
                # if measurement == 1:
                #     good_psi = prot.parity_projection_ket(e_psi, targets, 0, parity)
                #     # bad_psi = prot.parity_projection_ket(psi_initial, targets, 1, parity)
                # if measurement == -1:
                #     good_psi = prot.parity_projection_ket(e_psi, targets, 1, parity)
                #     # bad_psi = prot.parity_projection_ket(psi_initial, targets, 0, parity)


                # errs_permut += [prot.fidelity(rho, e_psi)]
                errs_permut += prot.fidelity(rho, e_psi)
            errs += [errs_permut]

        # Convert to arrays
        errs_array = np.array(errs)
        names_array = np.array(names)

        # Save arrays
        np.save("data/errors_" + subfix + err_type[0], errs_array)
        np.save("data/errors_names_" + subfix + err_type[0], names_array)


decompose(rho, good_psi, measurement, "good_")
decompose(rho, bad_psi, measurement, "bad_")


# Load all data
goodA = np.load("data/errors_good_A.npy")
badA = np.load("data/errors_bad_A.npy")
goodB = np.load("data/errors_good_B.npy")
badB = np.load("data/errors_bad_B.npy")
goodC = np.load("data/errors_good_C.npy")
badC = np.load("data/errors_bad_C.npy")
goodD = np.load("data/errors_good_D.npy")
badD = np.load("data/errors_bad_D.npy")


# The ones with D have a systematic err that needs to be removed
goodD = goodD[goodD != goodA.max()]
badD = badD[badD != badD.max()]

goods = 0
for i in [goodA, goodB, goodC, goodD]:
    goods += i.sum()

bads = 0
for i in [badA, badB, badC, badD]:
    bads += i.sum()
