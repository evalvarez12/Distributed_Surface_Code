"""
Decomposition routines into Kraus operators.

created-on: 30/06/17
"""
import protocols
import pauli_errors


# Set system size
system_size = 4

# Initialize protocols with errors
prot = protocols.Protocols(ps=0.0, pm=0.09, pg=0.09, pn=0.0 )

# The initial state
rho_initial = prot.operational_state(system_size)
# Parity measurement is made on one side of the bipartite state
targets = [0, 2, 4, 6]
measurement, rho = prot.monolithic(rho_initial, targets, "Z")

# Get all errors
pauli_errs = pauli_errors.get_pauli_errors(system_size, targets)

# for err_type in pauli_errs:
