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
prot = protocols.Protocols(ps=0.0, pm=0.0, pg=0.0, pn=0.0)

# The initial state
psi_initial = prot.operational_state_ket(system_size)
rho_initial = psi_initial * psi_initial.dag()

# Parity measurement is made on one side of the bipartite state
targets = [0, 2, 4, 6]
parity = "X"
measurement, rho = prot.expedient(rho_initial, targets, parity)
print("measurement: ", measurement)

# OLD for comparation
measurement_OLD, rho_OLD = prot.expedient_OLD(rho_initial, targets, parity)
print("measurement_OLD: ", measurement_OLD)



# Measurements to array to use .prod()
measurement = np.array(measurement)

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
def decompose(rho, psi, subfix):
    for err_type in pauli_errs:
        # Save the names of the errors
        names = err_type[1][1]

        # Array to store all fidelities results
        errs = []

        # Do the comparation with no pauli error
        if err_type[0] == "A":
            errs += [[prot.fidelity(rho, psi)]]
            names = ["I"] + names

        for permut in err_type[1][0]:
            errs_permut = []
            for e in permut:
                e_psi = e * psi
                errs_permut += [prot.fidelity(rho, e_psi)]
            errs += [errs_permut]

        # Convert to arrays
        errs_array = np.array([errs])
        names_array = np.array(names)

        # Save arrays
        np.save("data/errors_" + subfix + err_type[0], errs_array)
        np.save("data/errors_names_" + subfix + err_type[0], names_array)


decompose(rho, good_psi, "good_")
decompose(rho, bad_psi, "bad_")


good = np.load("data/errors_good_A.npy")
bad = np.load("data/errors_bad_A.npy")
