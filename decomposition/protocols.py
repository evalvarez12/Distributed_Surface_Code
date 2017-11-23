"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import circuit_block
import circuit
import error_models

# Determine parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.1
p_env = 2e-8
# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, pn, p_env)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

# First assemeble the small independent circuit
c_small = circuit.Circuit(p_env=p_env, circuit_block=cb.start_bell_pair)
c_small.add_circuit(circuit_block=cb.single_selection,
                    operation_qubits=[0, 1],
                    sigma="X")
c_small.add_circuit(circuit_block=cb.single_selection,
                    operation_qubits=[0, 1],
                    sigma="Z")

c_wrap_parallel = circuit.Circuit(p_env=p_env, circuit_block=c_small.run_parallel)




"""
Test Nickerson expedient protocol.
"""

print("------------------PROTOCOL 1-------------------")
c_pair_purification = circuit.Circuit(p_env=p_env, circuit_block=cb.start_bell_pair)
c_pair_purification.add_circuit(circuit_block=cb.double_selection,
                                operation_qubits=[0, 1],
                                sigma="X")
c_pair_purification.add_circuit(circuit_block=cb.double_selection,
                                operation_qubits=[0, 1],
                                sigma="Z")

c_ghz = circuit.Circuit(p_env=p_env,
                        circuit_block=c_pair_purification.run_parallel)
c_ghz.add_circuit(circuit_block=c_wrap_parallel.append_circuit)
c_ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                  targets=[1, 3, 0, 2], sigma="Z")
c_ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                  ancillas_pos=[4, 5, 6, 7],
                  projections=[0, 0, 0, 0])

# Get average number of steps
avg = 1
n_avg = 0
for i in range(avg):
    p, n, rho = c_ghz.run(None)
    n_avg += n
n_avg = n_avg/avg
print(rho)
print("n steps: ", n_avg)
print("F: ", qt.fidelity(rho, ghz_ref))


"""
Test Nickerson stringent protocol.
"""
print("------------------PROTOCOL 2-------------------")
# n, rho = cb.generate_bell_pair()
# # First two pumps of double selection
# cs = circuit.Circuit(p_env=p_env, circuit_block=cb.double_selection,
#                      operation_qubits=[0, 1],
#                      sigma="X")
# cs.add_circuit(circuit_block=cb.double_selection,
#                operation_qubits=[0, 1],
#                sigma="Z")
#
# # Append small circuit
# cs.add_circuit(circuit_block=c_small.run)
# # Two qubit gate
# cs.add_circuit(circuit_block=cb.two_qubit_gates,
#                controls=[2, 3], targets=[0, 1], sigma="Z")
# cs.add_circuit(circuit_block=cb.collapse_ancillas,
#                ancillas_pos=[2, 3], projections=[0, 0])
# # Append small circuit
# cs.add_circuit(circuit_block=c_small.run)
# # Two qubit gate
# cs.add_circuit(circuit_block=cb.two_qubit_gates,
#                controls=[2, 3], targets=[0, 1], sigma="X")# Get average number of steps
# cs.add_circuit(circuit_block=cb.collapse_ancillas,
#                ancillas_pos=[2, 3], projections=[0, 0])
#
# avg = 1
# n_avg = 0
# for i in range(avg):
#     p, n, rho = cs.run(rho)
#     n_avg += n
# n_avg = n_avg/avg
# print("n steps: ", n_avg)
# print("F: ", qt.fidelity(rho, rho_ref))
