import qutip as qt
import circuit_block
import circuit

# Determine parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.1
p_env = 5e-5
# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, pn, p_env)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()


print("------------------PROTOCOL 1-------------------")
_, _, rho = cb.start_bell_pair()
cs = circuit.Circuit(p_env=p_env, circuit_block=cb.single_selection,
                     operation_qubits=[0, 1],
                     sigma="X")
p, check, rho = cs.run(rho)
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------PROTOCOL 2-------------------")
_, _, rho = cb.start_bell_pair()
cs = circuit.Circuit(p_env=p_env, circuit_block=cb.single_selection,
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
_, _, rho = cb.start_bell_pair()
cs = circuit.Circuit(p_env=p_env, circuit_block=cb.double_selection,
                     operation_qubits=[0, 1],
                     sigma="X")

cs.add_circuit(circuit_block=cb.double_selection,
               operation_qubits=[0, 1],
               sigma="Z")
p, check, rho = cs.run(rho)
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))
