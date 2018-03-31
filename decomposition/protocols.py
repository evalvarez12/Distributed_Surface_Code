"""
Functions with the different protocols for generating a GHZ state.
They rely on circuit.py and circuit_block.py to enssemble the circuit that
generates the GHZ.
Each function returns a circuit object that must be executed to obtain a GHZ.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import circuit_block
import circuit


def pair_EPL(ps, pm, pg, eta, a0, a1, theta):
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_epl)
    return epl


def pair_single_sel(ps, pm, pg, eta, a0, a1, theta):
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)

    single_sel = circuit.Circuit(a0=a0, a1=a1,
                                 circuit_block=cb.single_selection_ops,
                                 targets=[0, 1], ancillas=[2, 3], sigma="Z")
    single_sel.add_circuit(circuit_block=epl.append_circuitMC)
    single_sel.add_circuit(circuit_block=cb.swap_pair,
                           pair=[0, 1])
    single_sel.add_circuit(circuit_block=cb.start_epl)

    return single_sel

def pair_double_sel(ps, pm, pg, eta, a0, a1, theta):
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)

    double_sel = circuit.Circuit(a0=a0, a1=a1,
                                 circuit_block=cb.double_selection_ops,
                                 targets=[0, 1], ancillas1=[2, 3],
                                 ancillas2=[4, 5], sigma="Z")
    double_sel.add_circuit(circuit_block=epl.append_circuitMC)
    double_sel.add_circuit(circuit_block=cb.swap_pair,
                           pair=[0, 1])
    double_sel.add_circuit(circuit_block=epl.append_circuitMC)
    double_sel.add_circuit(circuit_block=cb.swap_pair,
                           pair=[0, 1])
    double_sel.add_circuit(circuit_block=cb.start_epl)

    return double_sel






###################################################### OLD functions
#
# def EPL_4(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 4 created using 4 Bell pairs
#     generated using the EPL protocol.
#
#     Parameters
#     ----------
#     ps : (scalar) single qubit gate error rate.
#     pm : (scalar) measurement error rate.
#     pg : (scalar) two qubit gate error rate.
#     eta : (scalar) detection efficiency.
#     a0 : (scalar) extra environmental error when electron spin is being operated.
#     a1 : (scalar) default environmental error.
#     theta : (scalar) determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#     # First assemeble the small single selection circuit
#     add_EPL = circuit.Circuit(a0=a0, a1=a1,
#                               circuit_block=cb.start_epl)
#     wrap_EPL_parallel = circuit.Circuit(a0=a0, a1=a1,
#                                         circuit_block=add_EPL.run_parallel)
#
#     # Create the first initial EPL pair
#     start_EPL = circuit.Circuit(a0=a0, a1=a1,
#                                 circuit_block=cb.start_epl)
#
#     start_EPL.add_circuit(circuit_block=cb.swap_pair,
#                           pair=[0, 1])
#
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=start_EPL.run_parallel)
#
#     # Phase 2 - Create GHZ
#     # Create last two pairs
#     ghz.add_circuit(circuit_block=wrap_EPL_parallel.append_circuit)
#     # Apply two qubit gates in the nodes
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
#                     targets=[4, 5, 6, 7], sigma="X")
#     # Perform the measurements
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=4,
#                     measure_pos=[4, 5, 6, 7])
#     # Return the completed circuit
#     return ghz
#
# def EPL_4_simplified(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 4 created using 3 Bell pairs
#     generated using the EPL protocol.
#     Note: On simplified you need to substract the average time it
#     takes to create a EPL pair
#
#     Parameters
#     ----------
#     ps : (scalar) single qubit gate error rate.
#     pm : (scalar) measurement error rate.
#     pg : (scalar) two qubit gate error rate.
#     eta : (scalar) detection efficiency.
#     a0 : (scalar) extra environmental error when electron spin is being operated.
#     a1 : (scalar) default environmental error.
#     theta : (scalar) determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#     # First assemeble the small single selection circuit
#     start_EPL = circuit.Circuit(a0=a0, a1=a1,
#                                 circuit_block=cb.start_epl)
#
#     start_EPL.add_circuit(circuit_block=cb.swap_pair,
#                           pair=[0, 1])
#
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=start_EPL.run_parallel)
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=start_EPL.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
#                     targets=[4, 5], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=4,
#                     measure_pos=[4, 5])
#     # Return the completed circuit
#     return ghz
#
#
# def BK_4(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 4 created using 4 Bell pairs
#     generated using the BK protocol.
#
#     Parameters
#     ----------
#     ps : single qubit gate error rate.
#     pm : measurement error rate.
#     pg : two qubit gate error rate.
#     eta : detection efficiency.
#     a0 : extra environmental error when electron spin is being operated.
#     a1 : default environmental error.
#     theta : determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#
#     # First assemeble the small single selection circuit
#     start_BK = circuit.Circuit(a0=a0, a1=a1,
#                                circuit_block=cb.start_BK)
#     start_BK.add_circuit(circuit_block=cb.swap_pair,
#                          pair=[0, 1])
#
#     add_BK = circuit.Circuit(a0=a0, a1=a1,
#                              circuit_block=cb.start_BK)
#     # Wrapper circuit used to run in parallel the add_BK circuit
#     wrap_BK_parallel = circuit.Circuit(a0=a0, a1=a1,
#                                        circuit_block=add_BK.run_parallel)
#
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=start_BK.run_parallel)
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=wrap_BK_parallel.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
#                     targets=[4, 5, 6, 7], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=4,
#                     measure_pos=[4, 5, 6, 7])
#
#     # Return the completed circuit
#     return ghz
#
# def BK_4_simplified(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 4 created using 3 Bell pairs
#     generated using the BK protocol.
#     Note: On simplified you need to substract the average time it
#     takes to create a BK pair
#
#
#     Parameters
#     ----------
#     ps : single qubit gate error rate.
#     pm : measurement error rate.
#     pg : two qubit gate error rate.
#     eta : detection efficiency.
#     a0 : extra environmental error when electron spin is being operated.
#     a1 : default environmental error.
#     theta : determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#
#     start_BK = circuit.Circuit(a0=a0, a1=a1,
#                                circuit_block=cb.start_BK)
#     start_BK.add_circuit(circuit_block=cb.swap_pair,
#                          pair=[0, 1])
#     add_BK = circuit.Circuit(a0=a0, a1=a1,
#                              circuit_block=cb.start_BK)
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=start_BK.run_parallel)
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=add_BK.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
#                     targets=[4, 5], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=4,
#                     measure_pos=[4, 5])
#
#     # Return the completed circuit
#     return ghz
#
# def BK_3(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 3 created using 2 Bell pairs
#     generated using the BK protocol.
#
#     Parameters
#     ----------
#     ps : (scalar) single qubit gate error rate.
#     pm : (scalar) measurement error rate.
#     pg : (scalar) two qubit gate error rate.
#     eta : (scalar) detection efficiency.
#     a0 : (scalar) extra environmental error when electron spin is being operated.
#     a1 : (scalar) default environmental error.
#     theta : (scalar) determines how the states are initialized when generating remote
#             entanglement.
#     """
#
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#
#     # First assemeble the small single selection circuit
#     start_BK = circuit.Circuit(a0=a0, a1=a1,
#                                circuit_block=cb.start_BK)
#     start_BK.add_circuit(circuit_block=cb.swap_pair,
#                          pair=[0, 1])
#
#
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=cb.start_BK)
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=start_BK.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
#                     targets=[2], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=3,
#                     ancillas_pos=[2])
#
#     # Return the completed circuit
#     return ghz
#
# def EPL_3(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 3 created using 2 Bell pairs
#     generated using the EPL protocol.
#
#     Parameters
#     ----------
#     ps : (scalar) single qubit gate error rate.
#     pm : (scalar) measurement error rate.
#     pg : (scalar) two qubit gate error rate.
#     eta : (scalar) detection efficiency.
#     a0 : (scalar) extra environmental error when electron spin is being operated.
#     a1 : (scalar) default environmental error.
#     theta : (scalar) determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#
#     # First assemeble the small single selection circuit
#     start_EPL = circuit.Circuit(a0=a0, a1=a1,
#                                 circuit_block=cb.start_epl)
#     start_EPL.add_circuit(circuit_block=cb.swap_pair,
#                           pair=[0, 1])
#
#     # Phase 1 - Create initial pairs
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=cb.start_epl)
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=start_EPL.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
#                     targets=[2], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=3,
#                     ancillas_pos=[2])
#
#     # Return the completed circuit
#     return ghz
#
#
# def purification_simple_4(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 4 created using using the simple purification
#     protocol (see documentation).
#
#     Parameters
#     ----------
#     ps : (scalar) single qubit gate error rate.
#     pm : (scalar) measurement error rate.
#     pg : (scalar) two qubit gate error rate.
#     eta : (scalar) detection efficiency.
#     a0 : (scalar) extra environmental error when electron spin is being operated.
#     a1 : (scalar) default environmental error.
#     theta : (scalar) determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#
#     # First assemeble the small single selection circuit
#     single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
#                                          circuit_block=cb.start_epl)
#
#     single_sel_simple2.add_circuit(circuit_block=cb.single_selection,
#                                    operation_qubits=[0, 1],
#                                    sigma="Z")
#     # Wrapper used to run the circuit in parallel
#     wrap_single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
#                                               circuit_block=single_sel_simple2.run_parallel)
#
#     single_sel_simple1 = circuit.Circuit(a0=a0, a1=a1,
#                                          circuit_block=cb.start_epl)
#
#     single_sel_simple1.add_circuit(circuit_block=cb.single_selection,
#                                    operation_qubits=[0, 1],
#                                    sigma="Z")
#
#
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=single_sel_simple1.run_parallel)
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=wrap_single_sel_simple2.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
#                     targets=[4, 5, 6, 7], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=4,
#                     measure_pos=[4, 5, 6, 7])
#
#     # Return the completed circuit
#     return ghz
#
#
# def purification_simple_3(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 3 created using using the simple purification
#     protocol (see documentation).
#
#     Parameters
#     ----------
#     ps : (scalar) single qubit gate error rate.
#     pm : (scalar) measurement error rate.
#     pg : (scalar) two qubit gate error rate.
#     eta : (scalar) detection efficiency.
#     a0 : (scalar) extra environmental error when electron spin is being operated.
#     a1 : (scalar) default environmental error.
#     theta : (scalar) determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#
#     # First assemeble the small single selection circuit
#     single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
#                                          circuit_block=cb.start_epl)
#
#     single_sel_simple2.add_circuit(circuit_block=cb.single_selection,
#                                    operation_qubits=[0, 1],
#                                    sigma="Z")
#     # single_sel_simple2.add_circuit(circuit_block=cb.swap_pair,
#     #                                pair=[0, 1])
#
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1,
#                           circuit_block=cb.start_epl)
#     ghz.add_circuit(circuit_block=cb.single_selection,
#                     operation_qubits=[0, 1],
#                     sigma="Z")
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=single_sel_simple2.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
#                     targets=[2], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_Z,
#                     ancillas_pos=[2],
#                     projections=[0])
#
#     # Return the completed circuit
#     return ghz
#
#
# def purification_complex_4(ps, pm, pg, eta, a0, a1, theta):
#     """
#     GHZ state of weigth 4 created using using the complex purification
#     protocol (see documentation).
#
#     Parameters
#     ----------
#     ps : (scalar) single qubit gate error rate.
#     pm : (scalar) measurement error rate.
#     pg : (scalar) two qubit gate error rate.
#     eta : (scalar) detection efficiency.
#     a0 : (scalar) extra environmental error when electron spin is being operated.
#     a1 : (scalar) default environmental error.
#     theta : (scalar) determines how the states are initialized when generating remote
#             entanglement.
#     """
#     # Initialize the circuit block object
#     cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
#
#     # First assemeble the small independent circuit
#     double_sel1 = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_epl)
#     double_sel1.add_circuit(circuit_block=cb.swap_pair,
#                             pair=[0, 1])
#     double_sel1.add_circuit(circuit_block=cb.double_selection,
#                             operation_qubits=[0, 1],
#                             sigma="X")
#
#     double_sel2 = circuit.Circuit(a0=a0, a1=a1,
#                                   circuit_block=cb.start_epl)
#     double_sel2.add_circuit(circuit_block=cb.swap_pair,
#                             pair=[0, 1])
#     double_sel2.add_circuit(circuit_block=cb.double_selection,
#                             operation_qubits=[0, 1],
#                             sigma="X")
#     double_sel2.add_circuit(circuit_block=cb.swap_pair,
#                             pair=[0, 1])
#
#     wrap_double_sel2 = circuit.Circuit(a0=a0, a1=a1,
#                                        circuit_block=double_sel2.run_parallel)
#
#
#
#     # Phase 1 - Purify Bell pair
#     ghz = circuit.Circuit(a0=a0, a1=a1, circuit_block=double_sel1.run_parallel)
#
#
#     # Phase 2 - Create GHZ
#     ghz.add_circuit(circuit_block=wrap_double_sel2.append_circuit)
#     ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
#                     targets=[4, 5, 6, 7], sigma="X")
#     ghz.add_circuit(circuit_block=cb.collapse_ancillas_GHZ,
#                     ghz_size=4,
#                     measure_pos=[4, 5, 6, 7])
#
#     # Return the completed circuit
#     return ghz
