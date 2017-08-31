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
prot = protocols.Protocols(ps=0.0, pm=0.09, pg=0.09, pn=0.0)

# The initial state
state_initial = prot.operational_state_ket(system_size)
rho_initial = state_initial * state_initial.dag()

# Parity measurement is made on one side of the bipartite state
targets = [0, 2, 4, 6]
measurement, rho = prot.monolithic(rho_initial, targets, "Z")

# Measurements to array to use .prod()
measurement = np.array(measurement)

# Do the parity distintion to identify the good vs bad measurement
if measurement.prod() == 1:
    good_rho = projectors.project_even(rho)
    bad_rho = projectors.project_odd(rho)
if measurement.prod() == -1:
    good_rho = projectors.project_odd(rho)
    bad_rho = projectors.project_even(rho)


# Get all errors
pauli_errs = pauli_errors.get_pauli_errors(2*system_size, targets)


def decompose(rho, subfix):
    for err_type in pauli_errs:
        errs = []
        for permut in err_type[1][0]:
            errs_permut = []
            for e in permut:
                e_state = e * state_initial
                errs_permut += [prot.fidelity(rho, e_state)]
            errs += [errs_permut]

        # Convert to arrays
        errs_array = np.array([errs])
        names_array = np.array(err_type[1][1])

        # Save arrays
        np.save("data/errors_" + subfix + err_type[0], errs_array)
        np.save("data/errors_names_" + subfix + err_type[0], names_array)


print("measurement: ", measurement)
decompose(good_rho, "good_")
decompose(bad_rho, "bad_")
