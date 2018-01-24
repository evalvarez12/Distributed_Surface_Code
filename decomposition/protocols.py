"""
File to define the different protocols for generating a GHZ state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import circuit_block
import circuit


def EPL_4(ps, pm, pg, eta, a0, a1, theta):

    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    # First assemeble the small single selection circuit
    add_EPL = circuit.Circuit(a0=a0, a1=a1,
                              circuit_block=cb.start_epl)
    wrap_EPL_parallel = circuit.Circuit(a0=a0, a1=a1,
                                        circuit_block=add_EPL.run_parallel)

    start_EPL = circuit.Circuit(a0=a0, a1=a1,
                                circuit_block=cb.start_epl)

    start_EPL.add_circuit(circuit_block=cb.swap_pair,
                          pair=[0, 1])

    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=start_EPL.run_parallel)

    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=wrap_EPL_parallel.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                    targets=[1, 3, 0, 2], sigma="Z")
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
                    ancillas_pos=[4, 5, 6, 7],
                    projections=[0, 0, 0, 0])

    # Return the completed circuit
    return ghz

def EPL_4_simplified(ps, pm, pg, eta, a0, a1, theta):

    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    # First assemeble the small single selection circuit
    start_EPL = circuit.Circuit(a0=a0, a1=a1,
                                circuit_block=cb.start_epl)

    start_EPL.add_circuit(circuit_block=cb.swap_pair,
                          pair=[0, 1])

    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=start_EPL.run_parallel)

    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=start_EPL.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
                    targets=[4, 5], sigma="X")
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_Z,
                    ancillas_pos=[4, 5],
                    projections=[0, 0])

    # Return the completed circuit
    return ghz


def BK_4(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # First assemeble the small single selection circuit
    start_BK = circuit.Circuit(a0=a0, a1=a1,
                               circuit_block=cb.start_BK)
    start_BK.add_circuit(circuit_block=cb.swap_pair,
                         pair=[0, 1])

    add_BK = circuit.Circuit(a0=a0, a1=a1,
                             circuit_block=cb.start_BK)
    wrap_BK_parallel = circuit.Circuit(a0=a0, a1=a1,
                                       circuit_block=add_BK.run_parallel)

    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=start_BK.run_parallel)

    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=wrap_BK_parallel.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                    targets=[1, 3, 0, 2], sigma="Z")
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
                    ancillas_pos=[4, 5, 6, 7],
                    projections=[0, 0, 0, 0])

    # Return the completed circuit
    return ghz

def BK_4_simplified(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    start_BK = circuit.Circuit(a0=a0, a1=a1,
                               circuit_block=cb.start_BK)
    start_BK.add_circuit(circuit_block=cb.swap_pair,
                         pair=[0, 1])

    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=start_BK.run_parallel)

    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=start_BK.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
                    targets=[4, 5], sigma="X")
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_Z,
                    ancillas_pos=[4, 5],
                    projections=[0, 0])

    # Return the completed circuit
    return ghz

def BK_3(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # First assemeble the small single selection circuit
    start_BK = circuit.Circuit(a0=a0, a1=a1,
                               circuit_block=cb.start_BK)
    start_BK.add_circuit(circuit_block=cb.swap_pair,
                         pair=[0, 1])


    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_BK)

    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=start_BK.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
                    targets=[2], sigma="X")
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_Z,
                    ancillas_pos=[2],
                    projections=[0])

    # Return the completed circuit
    return ghz

def EPL_3(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # First assemeble the small single selection circuit
    start_EPL = circuit.Circuit(a0=a0, a1=a1,
                               circuit_block=cb.start_epl)
    start_EPL.add_circuit(circuit_block=cb.swap_pair,
                         pair=[0, 1])


    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_epl)

    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=start_EPL.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
                    targets=[2], sigma="X")
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_Z,
                    ancillas_pos=[2],
                    projections=[0])

    # Return the completed circuit
    return ghz


def purification_simple_4(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # First assemeble the small single selection circuit
    single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=cb.start_epl)

    single_sel_simple2.add_circuit(circuit_block=cb.single_selection,
                                   operation_qubits=[0, 1],
                                   sigma="X")
    # single_sel_simple2.add_circuit(circuit_block=cb.swap_pair,
    #                                pair=[0, 1])
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
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
                    ancillas_pos=[4, 5, 6, 7],
                    projections=[0, 0, 0, 0])

    # Return the completed circuit
    return ghz


def purification_simple_3(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # First assemeble the small single selection circuit
    single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=cb.start_epl)

    single_sel_simple2.add_circuit(circuit_block=cb.single_selection,
                                   operation_qubits=[0, 1],
                                   sigma="X")
    # single_sel_simple2.add_circuit(circuit_block=cb.swap_pair,
    #                                pair=[0, 1])

    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_epl)
    ghz.add_circuit(circuit_block=cb.single_selection,
                    operation_qubits=[0, 1],
                    sigma="X")

    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=single_sel_simple2.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
                    targets=[2], sigma="X")
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_Z,
                    ancillas_pos=[2],
                    projections=[0])

    # Return the completed circuit
    return ghz


def purification_medium_4(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # First assemeble the small independent circuit
    single_sel_medium2 = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=cb.start_epl)
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
                                         circuit_block=cb.start_epl)
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
    ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
                    ancillas_pos=[4, 5, 6, 7],
                    projections=[0, 0, 0, 0])

    # Return the completed circuit
    return ghz

def purification_complex_4(ps, pm, pg, eta, a0, a1, theta):
    # Initialize the circuit block object
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # First assemeble the small independent circuit
    double_sel1 = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_epl)
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
    # double_sel1.add_circuit(circuit_block=cb.collapse_ancillas_X,
    #                         ancillas_pos=[2, 3],
    #                         projections=[0, 0])
    # double_sel1.add_circuit(circuit_block=single_sel_medium2.append_circuit)
    # double_sel1.add_circuit(circuit_block=cb.two_qubit_gates, controls=[2, 3],
    #                         targets=[0, 1], sigma="X")
    # double_sel1.add_circuit(circuit_block=cb.collapse_ancillas_X,
    #                         ancillas_pos=[2, 3],
    #                         projections=[0, 0])

    single_sel_medium1 = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=cb.start_epl)
    single_sel_medium1.add_circuit(circuit_block=cb.swap_pair,
                                   pair=[0, 1])
    single_sel_medium1.add_circuit(circuit_block=cb.single_selection,
                                   operation_qubits=[0, 1],
                                   sigma="Z")
    single_sel_medium1.add_circuit(circuit_block=cb.single_selection,
                                   operation_qubits=[0, 1],
                                   sigma="X")

    single_sel_medium2 = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=cb.start_epl)
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



    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1, circuit_block=double_sel1.run_parallel)


    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=wrap_single_sel_medium2.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                    targets=[1, 3, 0, 2], sigma="Z")

    ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
                    ancillas_pos=[4, 5, 6, 7],
                    projections=[0, 0, 0, 0])

    # Return the completed circuit
    return ghz
