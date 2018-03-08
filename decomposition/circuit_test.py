"""
Test routines to test the functions in 'circuit.py'

author: Eduardo Villasenor
created-on: 20/11/17
"""
import qutip as qt
import circuit_block
import circuit
import numpy as np

# Determine parameters
ps = 0.005
pm = 0.005
pg = 0.005
a0 = 0
a1 = 0.
eta = 1/10.
theta = np.pi/4

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
cb_ideal = circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4)

print("------------------PROTOCOL 1-------------------")
_, _, rho = cb.start_epl()
print("Initial F", qt.fidelity(rho, qt.bell_state('10')))

cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb_ideal.single_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
rho, check = cs.run(rho)
print("check: ", check)
print("F: ", qt.fidelity(rho, qt.bell_state('10')))

# print("------------------PROTOCOL 2-------------------")
# _, _, rho = cb.start_epl()
# cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.single_selection,
#                      operation_qubits=[0, 1],
#                      sigma="X")
# cs.add_circuit(circuit_block=cb.single_selection,
#                operation_qubits=[0, 1],
#                sigma="Z")
# rho, check = cs.run(rho)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, qt.bell_state('10')))
#
#
# print("------------------PROTOCOL 3-------------------")
# _, _, rho = cb.start_epl()
# cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.double_selection,
#                      operation_qubits=[0, 1],
#                      sigma="X")
#
# cs.add_circuit(circuit_block=cb.double_selection,
#                operation_qubits=[0, 1],
#                sigma="Z")
# rho, check = cs.run(rho)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, qt.bell_state('10')))
