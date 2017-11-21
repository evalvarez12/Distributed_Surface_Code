import qutip as qt
import circuit_block
import circuit
import error_models as errs
import operations as ops
import numpy as np

# Determine parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.1

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, pn)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()


print("------------------PROTOCOL 1-------------------")
n, rho = cb.generate_bell_pair()
cs = circuit.Circuit(cb.single_selection, operation_qubits=[0, 1],
                     sigma="X")
n, rho = cs.run(rho)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------PROTOCOL 2-------------------")
n, rho = cb.generate_bell_pair()
cs = circuit.Circuit(circuit_block=cb.single_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
cs.add_circuit(link=True, circuit_block=cb.single_selection,
               operation_qubits=[0, 1],
               sigma="Z")
n, rho = cs.run(rho)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))


print("------------------PROTOCOL 3-------------------")
n, rho = cb.generate_bell_pair()
cs = circuit.Circuit(circuit_block=cb.double_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
cs.add_circuit(link=True, circuit_block=cb.double_selection,
               operation_qubits=[0, 1],
               sigma="Z")
n, rho = cs.run(rho)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))


print("------------------PROTOCOL 4-------------------")
n, rho = cb.generate_bell_pair()
cs = circuit.Circuit(circuit_block=cb.double_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
# NOTE: Here link=False makes no sense in a real circuit,
# as failure of the purification destroys the state rho
cs.add_circuit(link=False, circuit_block=cb.double_selection,
               operation_qubits=[0, 1],
               sigma="Z")
n, rho = cs.run(rho)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))
