import qutip as qt
import circuits_blocks
import circuits
import error_models as errs
import operations as ops
import numpy as np

# Determine parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.1

# Initialize  objects
cb = circuits_blocks.Blocks(ps, pm, pg, pn)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()


print("------------------PROTOCOL 1-------------------")
n, rho = cb.generate_bell_pair()
cs = circuits.Circuit(cb.single_selection, operation_qubits=[0, 1],
                                            sigma="X")
n, rho = cs.run(rho, 0, 0)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------PROTOCOL 2-------------------")
n, rho = cb.generate_bell_pair()
cs = circuits.Circuit(circuit_block=cb.single_selection,
                      operation_qubits=[0, 1],
                      sigma="X")
cs.add_circuit(dependant=True, circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="Z")
n, rho = cs.run(rho, 0, 0)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))




#
# print("------------------SINGLE SELECTION-------------------")
# n, rho_i = cb.generate_bell_pair()
# p, n, rho = cb.single_selection(rho_i, [0, 1], "X")
# print("p_success: ", p)
# print("n steps: ", n)
# print("F: ", qt.fidelity(rho, rho_i))
#
# print("------------------DOUBLE SELECTION-------------------")
# n, rho_i = cb.generate_bell_pair()
# p, n, rho = cb.double_selection(rho_i, [0, 1], "X")
# print("p_success: ", p)
# print("n steps: ", n)
# print("F: ", qt.fidelity(rho, rho_i))
