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
ps = 0.009
pm = 0.009
pg = 0.009
a0 = 3.
a1 = 1/80.
eta = 1/100.
theta = .24

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
cb_ideal = circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4)

print("------------------PROTOCOL 1-------------------")
epl = circuit.Circuit(a0=a0, a1=a1,
                      circuit_block=cb.start_epl)
rho, c = epl.run()
print("Initial F", qt.fidelity(rho, qt.bell_state('00')))
print(c)


print("------------------PROTOCOL 2-------------------")
epl = circuit.Circuit(a0=a0, a1=a1,
                      circuit_block=cb.start_epl)

circ = circuit.Circuit(a0=a0, a1=a1,
                       circuit_block=epl.run_parallel)
rho, c = circ.run()
print(rho)
print(c)

# _, _, rho = cb.start_epl()
# H = qt.snot(2, 0) * qt.snot(2, 1)
# # rho = H * rho * H.dag()
# cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.single_selection,
#                      operation_qubits=[0, 1],
#                      sigma="Z")
# # cs.add_circuit(circuit_block=cb.single_selection,
# #                operation_qubits=[0, 1],
# #                sigma="X")
# rho, check = cs.run(rho)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, qt.bell_state('00')))
#
#
# print("------------------PROTOCOL 3-------------------")
# _, _, rho = cb.start_epl()
# H = qt.snot(2, 0) * qt.snot(2, 1)
# # rho = H * rho * H.dag()
# cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.double_selection,
#                      operation_qubits=[0, 1],
#                      sigma="X")
#
# # cs.add_circuit(circuit_block=cb.double_selection,
# #                operation_qubits=[0, 1],
# #                sigma="X")
# rho, check = cs.run(rho)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, qt.bell_state('00')))
