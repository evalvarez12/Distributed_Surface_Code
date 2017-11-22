"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import circuit_block
import circuit

# Determine parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.1

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, pn)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()



"""
Test Nickerson expedient protocol.
"""
print("------------------PROTOCOL 1-------------------")
n, rho = cb.generate_bell_pair()
cs = circuit.Circuit(circuit_block=cb.double_selection,
                      operation_qubits=[0, 1],
                      sigma="X")
cs.add_circuit(circuit_block=cb.double_selection,
                               operation_qubits=[0, 1],
                               sigma="Z")
# Get average number of steps
avg = 1
n_avg = 0
for i in range(avg):
    p, n, rho = cs.run(rho)
    n_avg += n
n_avg = n_avg/avg
print("n steps: ", n_avg)
print("F: ", qt.fidelity(rho, rho_ref))


"""
Test Nickerson stringent protocol.
"""
print("------------------PROTOCOL 2-------------------")
# First assemeble the small independent circuit
c_small = circuit.Circuit(circuit_block=cb.add_bell_pair)
c_small.add_circuit(circuit_block=cb.single_selection,
                    operation_qubits=[2, 3],
                    sigma="X")
c_small.add_circuit(circuit_block=cb.single_selection,
                    operation_qubits=[2, 3],
                    sigma="Z")



n, rho = cb.generate_bell_pair()
# First two pumps of double selection
cs = circuit.Circuit(circuit_block=cb.double_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
cs.add_circuit(circuit_block=cb.double_selection,
               operation_qubits=[0, 1],
               sigma="Z")

# Append small circuit
cs.add_circuit(circuit_block=c_small.run)
# Two qubit gate
cs.add_circuit(circuit_block=cb.two_qubit_gates,
               controls=[2, 3], targets=[0, 1], sigma="Z")
cs.add_circuit(circuit_block=cb.collapse_ancillas,
               ancillas_pos=[2, 3], projections=[0, 0])
# Append small circuit
cs.add_circuit(circuit_block=c_small.run)
# Two qubit gate
cs.add_circuit(circuit_block=cb.two_qubit_gates,
               controls=[2, 3], targets=[0, 1], sigma="X")# Get average number of steps
cs.add_circuit(circuit_block=cb.collapse_ancillas,
               ancillas_pos=[2, 3], projections=[0, 0])

avg = 1
n_avg = 0
for i in range(avg):
    p, n, rho = cs.run(rho)
    n_avg += n
n_avg = n_avg/avg
print("n steps: ", n_avg)
print("F: ", qt.fidelity(rho, rho_ref))
