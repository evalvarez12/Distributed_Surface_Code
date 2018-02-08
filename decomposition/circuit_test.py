"""
Test routines to test the functions in 'circuit.py'

author: Eduardo Villasenor
created-on: 20/11/17
"""
import qutip as qt
import circuit_block
import circuit

# Determine parameters
ps = 0.006
pm = 0.006
pg = 0.006
a0 = 0
a1 = 0.
eta = 1
theta = .24

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()


print("------------------PROTOCOL 1-------------------")
_, _, rho = cb.start_epl()
cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.single_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
p, check, rho = cs.run(rho)
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------PROTOCOL 2-------------------")
_, _, rho = cb.start_epl()
cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.single_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
cs.add_circuit(circuit_block=cb.single_selection,
               operation_qubits=[0, 1],
               sigma="Z")
p, check, rho = cs.run(rho)
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))


print("------------------PROTOCOL 3-------------------")
_, _, rho = cb.start_epl()
cs = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.double_selection,
                     operation_qubits=[0, 1],
                     sigma="X")

cs.add_circuit(circuit_block=cb.double_selection,
               operation_qubits=[0, 1],
               sigma="Z")
p, check, rho = cs.run(rho)
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))
