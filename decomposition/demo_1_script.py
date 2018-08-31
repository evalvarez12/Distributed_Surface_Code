"""
    This is just Demo 1 as a script, if there are differences between
    the notebook and this, then trust the notebook.
"""

import qutip as qt
import numpy as np
import circuit
import circuit_block
import stabilizer

# Set parameters
# Gate and measurement error rates
ps = 0.003
pm = 0.003
pg = 0.003

# Very optimistic enviromental error rate and entanglement generation rate
a0 = 1.5
a1 = 1 / 80.
eta = 1 / 100.

# Theta to optimize entanglement generation
theta = .24

# Initialize the circuit block and stabilizer objects 
# Circuits are assemebled in reversed order due to recursion
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
epl = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_epl)

# Phase 3 - Create GHZ
# Perform the measurements
protocol = circuit.Circuit(a0=a0, a1=a1,
                           circuit_block=cb.collapse_ancillas_GHZ,
                           ghz_size=4,
                           measure_pos=[4, 5, 6, 7])
# Apply two qubit gates in the nodes
protocol.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
                     targets=[4, 5, 6, 7], sigma="X")

# Phase 2 Create last two pairs
pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                       pair=[0, 1])
pair.add_circuit(circuit_block=cb.start_epl)
# Wrapper to run in parallel
wrap_EPL_parallel = circuit.Circuit(a0=a0, a1=a1,
                                    circuit_block=pair.run_parallel)

protocol.add_circuit(circuit_block=wrap_EPL_parallel.append_circuit)

# Phase 1 Create initial two pairs

protocol.add_circuit(circuit_block=epl.run_parallel)

ghz, resources = protocol.run()

# Reference GHZ state to compare the result
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

print("Resulting state fidelity: ", qt.fidelity(ghz, ghz_ref))
print("Resources used: ", resources)

# Initialize stabilizer object
stab = stabilizer.Stabilizer(ps=ps, pm=pm, pg=pg)

# Lists to save fidelities and times
fidelity = []
times = []

iterations = 5
rho = ghz * 0
for i in range(iterations):
    r, c = protocol.run()

    # Twirl the resulting state to distribute any error symmetrically over the qubits
    r = stab.twirl_ghz(r)
    
    fidelity += [qt.fidelity(r, ghz_ref)]
    rho += r
    times += [c["time"]]
    
rho = rho/iterations
print("Fidelity: ", np.average(fidelity), np.std(fidelity))
print("Time:", np.average(times), np.std(times))

import noise_modeling

# Select a type of stabilizer to extract the model from
stab_type = "X"
stab_size = 4

# Initialize objects
model = noise_modeling.NoiseModel(stab_size, stab_type)
model.separate_basis_parity()

# Choi state for noise noise modeling
choi = model._choi_state_ket(stab_size)
choi = choi * choi.dag()
targets = list(range(stab_size))

probs, rhos = stab.measure_ghz_stabilizer(choi, ghz, targets, stab_type)

# Set channel output and make chi matrix
model.set_rho(rhos, probs)
model.make_chi_matrix()
    
# Sanity check: Check the total sum of the decomposition = 1
print("Sum of all probabilties: ", model.check_total_sum())

print("Decomposition: ", model.chi)

I_OK = model.chi["IIII_OK"]
I_NOK = model.chi["IIII_NOK"]
# The sum of all physical errors
E = (1 - model.chi["IIII_OK"] - model.chi["IIII_NOK"])/4.
print("Errors:")
print("OK: ", I_OK)
print("NOK: ", I_NOK)
print("QUBIT ERROR: ", E)

import tools.names as names
import pickle


file_name = names.chi(ps, pm, pg, eta, a0, a1, theta,
                      stab_size, stab_type, "DEMO")

print("File name: ", file_name)
pickle_out = open(file_name, "wb")
pickle.dump(model.chi, pickle_out, protocol=2)
pickle_out.close()