"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import collections
import circuit_block
import circuit
import error_models

# Determine parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 83.33
a1 = 1/3.
eta = (0.1)*(0.03)*(0.8)
theta = np.pi/4.

iterations = 200

# Initialize objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

# Lists to save results
SIMPLE_fidelity = []
SIMPLE_steps = []
MEDIUM_fidelity = []
MEDIUM_steps = []
COMPLEX_fidelity = []
COMPLEX_steps = []

'''PROTOCOLS START HERE'''

print("-------------------PROTOCOL SIMPLE------------------")
# First assemeble the small single selection circuit
single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
                                     circuit_block=cb.start_epl)

single_sel_simple2.add_circuit(circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="X")
single_sel_simple2.add_circuit(circuit_block=cb.swap_pair,
                               pair=[0, 1])
wrap_single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
                                          circuit_block=single_sel_simple2.run_parallel)


single_sel_simple1 = circuit.Circuit(a0=a0, a1=a1,
                                     circuit_block=cb.start_epl)

single_sel_simple1.add_circuit(circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="X")


# Phase 1 - Purify Bell pair
ghz = circuit.Circuit(a0=a0, a1=a1,
                      circuit_block=single_sel_simple1.run_parallel)

# Phase 2 - Create GHZ
ghz.add_circuit(circuit_block=wrap_single_sel_simple2.append_circuit)
ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                targets=[1, 3, 0, 2], sigma="Z")
ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                ancillas_pos=[4, 5, 6, 7],
                projections=[0, 0, 0, 0])

# Get average number of steps
fidelity = []
check = collections.Counter({})
for i in range(iterations):
    p, c, rho = ghz.run(None)
    check += c
    fidelity += [qt.fidelity(rho, ghz_ref)]

# print(rho)
SIMPLE_fidelity += [(np.average(fidelity), np.std(fidelity))]
# SIMPLE_steps += [(np.average(steps), np.std(steps))]
print("End protocol")

# Average over check values
for k in check:
    check[k] = check[k]/iterations
print("check: ", check)
print("F: ", np.average(fidelity), np.std(fidelity))



"""
Nickerson expedient protocol.
"""
print("------------------PROTOCOL MEDIUM-------------------")
# First assemeble the small independent circuit
single_sel_medium2 = circuit.Circuit(a0=a0, a1=a1,
                                     circuit_block=cb.start_bell_pair)
single_sel_medium2.add_circuit(circuit_block=cb.swap_pair,
                               pair=[0, 1])
single_sel_medium2.add_circuit(circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="Z")
single_sel_medium2.add_circuit(circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="X")
single_sel_medium2.add_circuit(circuit_block=cb.swap_pair,
                               pair=[0, 1])

wrap_single_sel_medium2 = circuit.Circuit(a0=a0, a1=a1,
                                          circuit_block=single_sel_medium2.run_parallel)


single_sel_medium1 = circuit.Circuit(a0=a0, a1=a1,
                                     circuit_block=cb.start_bell_pair)
single_sel_medium1.add_circuit(circuit_block=cb.swap_pair,
                               pair=[0, 1])
single_sel_medium1.add_circuit(circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="Z")
single_sel_medium1.add_circuit(circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="X")

# Phase 1 - Purify Bell pair
ghz = circuit.Circuit(a0=a0, a1=a1,
                      circuit_block=single_sel_medium1.run_parallel)

# Phase 2 - Create GHZ
ghz.add_circuit(circuit_block=wrap_single_sel_medium2.append_circuit)
ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                targets=[1, 3, 0, 2], sigma="Z")
ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                ancillas_pos=[4, 5, 6, 7],
                projections=[0, 0, 0, 0])

# Get average number of steps
fidelity = []
check = collections.Counter({})
for i in range(iterations):
    p, c, rho = ghz.run(None)
    check += c
    fidelity += [qt.fidelity(rho, ghz_ref)]

for k in check:
    check[k] = check[k]/iterations

# print(rho)
MEDIUM_fidelity += [(np.average(fidelity), np.std(fidelity))]
# MEDIUM_steps += [(np.average(steps), np.std(steps))]
print("End protocol")


print("check: ", check)
print("F: ", np.average(fidelity), np.std(fidelity))


"""
Nickerson stringent protocol.
"""
print("------------------PROTOCOL COMPLEX-------------------")
# First assemeble the small independent circuit
double_sel1 = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_bell_pair)
double_sel1.add_circuit(circuit_block=cb.swap_pair,
                        pair=[0, 1])
double_sel1.add_circuit(circuit_block=cb.double_selection,
                        operation_qubits=[0, 1],
                        sigma="Z")
double_sel1.add_circuit(circuit_block=cb.double_selection,
                        operation_qubits=[0, 1],
                        sigma="X")
# double_sel1.add_circuit(circuit_block=single_sel_medium2.append_circuit)
# double_sel1.add_circuit(circuit_block=cb.two_qubit_gates, controls=[2, 3],
#                         targets=[0, 1], sigma="Z")
# double_sel1.add_circuit(circuit_block=cb.collapse_ancillas,
#                         ancillas_pos=[2, 3],
#                         projections=[0, 0])
# double_sel1.add_circuit(circuit_block=single_sel_medium2.append_circuit)
# double_sel1.add_circuit(circuit_block=cb.two_qubit_gates, controls=[2, 3],
#                         targets=[0, 1], sigma="X")
# double_sel1.add_circuit(circuit_block=cb.collapse_ancillas,
#                         ancillas_pos=[2, 3],
#                         projections=[0, 0])

# Phase 1 - Purify Bell pair
ghz = circuit.Circuit(a0=a0, a1=a1, circuit_block=double_sel1.run_parallel)


# Phase 2 - Create GHZ
ghz.add_circuit(circuit_block=wrap_single_sel_medium2.append_circuit)
ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                targets=[1, 3, 0, 2], sigma="Z")

ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                ancillas_pos=[4, 5, 6, 7],
                projections=[0, 0, 0, 0])
# ghz.add_circuit(circuit_block=wrap_single_sel_medium2.append_circuit)
# ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
#                 targets=[1, 3, 0, 2], sigma="Z")
#
# ghz.add_circuit(circuit_block=cb.collapse_ancillas,
#                 ancillas_pos=[4, 5, 6, 7],
#                 projections=[0, 0, 0, 0])

# Get average number of steps
fidelity = []
check = collections.Counter({})
for i in range(iterations):
    p, c, rho = ghz.run(None)
    check += c
    fidelity += [qt.fidelity(rho, ghz_ref)]

COMPLEX_fidelity += [(np.average(fidelity), np.std(fidelity))]

for k in check:
    check[k] = check[k]/iterations

# COMPLEX_steps += [(np.average(steps), np.std(steps))]
# print(rho)
print("End protocol")
print("check: ", check)
print("F: ", np.average(fidelity), np.std(fidelity))

# np.save("data/SIMPLE_fidelity_BK", SIMPLE_fidelity)
# np.save("data/SIMPLE_steps_BK", SIMPLE_steps)
# np.save("data/MEDIUM_fidelity_BK", MEDIUM_fidelity)
# np.save("data/MEDIUM_steps_BK", MEDIUM_steps)
# np.save("data/COMPLEX_fidelity_BK", COMPLEX_fidelity)
# np.save("data/COMPLEX_steps_BK", COMPLEX_steps)
# np.save("data/p_env", p_env_var)
