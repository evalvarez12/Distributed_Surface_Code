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
p_env = 5e-4

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, pn, p_env)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

# Lists to save results
SIMPLE_fidelity = []
SIMPLE_steps = []
MEDIUM_fidelity = []
MEDIUM_steps = []
COMPLEX_fidelity = []
COMPLEX_steps = []

p_env_var = np.geomspace(5e-3, 5e-6, 20)

for p_env in p_env_var:
    cb.change_parameters(ps, pm, pg, pn, p_env)

    # print("-------------------PROTOCOL SIMPLE------------------")
    # First assemeble the small single selection circuit
    c_single_sel1 = circuit.Circuit(p_env=p_env, circuit_block=cb.start_bell_pair)
    c_single_sel1.add_circuit(circuit_block=cb.single_selection,
                              operation_qubits=[0, 1],
                              sigma="X")
    c_wrap_single_sel1 = circuit.Circuit(p_env=p_env,
                                         circuit_block=c_single_sel1.run_parallel)

    c_ghz = circuit.Circuit(p_env=p_env,
                            circuit_block=c_single_sel1.run_parallel)
    c_ghz.add_circuit(circuit_block=c_wrap_single_sel1.append_circuit)

    # Phase 2 - Create GHZ
    c_ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                      targets=[1, 3, 0, 2], sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                      ancillas_pos=[4, 5, 6, 7],
                      projections=[0, 0, 0, 0])

    # Get average number of steps
    avg = 1000
    fidelity = []
    steps = []
    for i in range(avg):
        p, n, rho = c_ghz.run(None)
        steps += [n]
        fidelity += [qt.fidelity(rho, ghz_ref)]

    # print(rho)
    SIMPLE_fidelity += [(np.average(fidelity), np.std(fidelity))]
    SIMPLE_steps += [(np.average(steps), np.std(steps))]
    # print("End protocol")
    # print("n steps: ", np.average(steps), np.std(steps))
    # print("F: ", np.average(fidelity), np.std(fidelity))



    """
    Nickerson expedient protocol.
    """
    # print("------------------PROTOCOL MEDIUM-------------------")
    # First assemeble the small independent circuit
    c_single_sel2 = circuit.Circuit(p_env=p_env, circuit_block=cb.start_bell_pair)
    c_single_sel2.add_circuit(circuit_block=cb.single_selection,
                              operation_qubits=[0, 1],
                              sigma="X")
    c_single_sel2.add_circuit(circuit_block=cb.single_selection,
                              operation_qubits=[0, 1],
                              sigma="Z")

    c_wrap_single_sel2 = circuit.Circuit(p_env=p_env,
                                         circuit_block=c_single_sel2.run_parallel)

    c_double_sel = circuit.Circuit(p_env=p_env,
                                   circuit_block=cb.start_bell_pair)
    c_double_sel.add_circuit(circuit_block=cb.single_selection,
                             operation_qubits=[0, 1],
                             sigma="X")
    c_double_sel.add_circuit(circuit_block=cb.single_selection,
                             operation_qubits=[0, 1],
                             sigma="Z")

    c_ghz = circuit.Circuit(p_env=p_env,
                            circuit_block=c_double_sel.run_parallel)

    # Phase 2 - Create GHZ
    c_ghz.add_circuit(circuit_block=c_wrap_single_sel2.append_circuit)
    c_ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                      targets=[1, 3, 0, 2], sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                      ancillas_pos=[4, 5, 6, 7],
                      projections=[0, 0, 0, 0])

    # Get average number of steps
    avg = 1000
    fidelity = []
    steps = []
    for i in range(avg):
        p, n, rho = c_ghz.run(None)
        steps += [n]
        fidelity += [qt.fidelity(rho, ghz_ref)]

    # print(rho)
    MEDIUM_fidelity += [(np.average(fidelity), np.std(fidelity))]
    MEDIUM_steps += [(np.average(steps), np.std(steps))]
    # print("End protocol")
    # print("n steps: ", np.average(steps), np.std(steps))
    # print("F: ", np.average(fidelity), np.std(fidelity))


    """
    Nickerson stringent protocol.
    """
    # print("------------------PROTOCOL COMPLEX-------------------")
    # Reuse pair purification circuit from previous protocol
    # and single selection parallel
    c_ghz = circuit.Circuit(p_env=p_env,
                            circuit_block=c_double_sel.run_parallel)
    c_ghz.add_circuit(circuit_block=c_wrap_single_sel2.append_circuit)
    c_ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                      targets=[0, 1, 2, 3], sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.single_selection, operation_qubits=[4, 5],
                      sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.single_selection, operation_qubits=[6, 7],
                      sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                      ancillas_pos=[4, 5, 6, 7],
                      projections=[0, 0, 0, 0])

    # Phase 2 - Create GHZ
    c_ghz.add_circuit(circuit_block=c_wrap_single_sel2.append_circuit)
    c_ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                      targets=[1, 3, 0, 2], sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.single_selection, operation_qubits=[4, 5],
                      sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.single_selection, operation_qubits=[6, 7],
                      sigma="Z")
    c_ghz.add_circuit(circuit_block=cb.collapse_ancillas,
                      ancillas_pos=[4, 5, 6, 7],
                      projections=[0, 0, 0, 0])

    # Get average number of steps
    avg = 1000
    fidelity = []
    steps = []
    for i in range(avg):
        p, n, rho = c_ghz.run(None)
        steps += [n]
        fidelity += [qt.fidelity(rho, ghz_ref)]

    COMPLEX_fidelity += [(np.average(fidelity), np.std(fidelity))]
    COMPLEX_steps += [(np.average(steps), np.std(steps))]
    # print(rho)
    # print("End protocol")
    # print("n steps: ", np.average(steps), np.std(steps))
    # print("F: ", np.average(fidelity), np.std(fidelity))

np.save("data/SIMPLE_fidelity_BK", SIMPLE_fidelity)
np.save("data/SIMPLE_steps_BK", SIMPLE_steps)
np.save("data/MEDIUM_fidelity_BK", MEDIUM_fidelity)
np.save("data/MEDIUM_steps_BK", MEDIUM_steps)
np.save("data/COMPLEX_fidelity_BK", COMPLEX_fidelity)
np.save("data/COMPLEX_steps_BK", COMPLEX_steps)
np.save("data/p_env", p_env_var)
