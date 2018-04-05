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
    """
    Create a entangled pair using the EPL protocol.

    Parameters
    ----------
    ps : single qubit gate error rate.
    pm : measurement error rate.
    pg : two qubit gate error rate.
    eta : detection efficiency.
    a0 : extra environmental error when electron spin is being operated.
    a1 : default environmental error.
    theta : determines how the states are initialized when generating remote
            entanglement.
    """
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_epl)
    return epl


def pair_single_sel(ps, pm, pg, eta, a0, a1, theta):
    """
    Create a entangled pair using single selection
    distillation protocol. The EPL protocol is used in the intermediate
    entangled pairs used.

    Parameters
    ----------
    ps : single qubit gate error rate.
    pm : measurement error rate.
    pg : two qubit gate error rate.
    eta : detection efficiency.
    a0 : extra environmental error when electron spin is being operated.
    a1 : default environmental error.
    theta : determines how the states are initialized when generating remote
            entanglement.
    """
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)

    single_sel = circuit.Circuit(a0=a0, a1=a1,
                                 circuit_block=cb.single_selection_ops,
                                 targets=[0, 1], ancillas=[2, 3], sigma="Z")
    single_sel.add_circuit(circuit_block=cb.swap_pair,
                           pair=[0, 1])
    single_sel.add_circuit(circuit_block=epl.append_circuit)
    single_sel.add_circuit(circuit_block=cb.start_epl)

    return single_sel

def pair_double_sel(ps, pm, pg, eta, a0, a1, theta):
    """
    Create a entangled pair using the double selection
    distillation protocol. The EPL protocol is used in the intermediate
    entangled pairs used.

    Parameters
    ----------
    ps : single qubit gate error rate.
    pm : measurement error rate.
    pg : two qubit gate error rate.
    eta : detection efficiency.
    a0 : extra environmental error when electron spin is being operated.
    a1 : default environmental error.
    theta : determines how the states are initialized when generating remote
            entanglement.
    """
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)

    double_sel = circuit.Circuit(a0=a0, a1=a1,
                                 circuit_block=cb.double_selection_ops,
                                 targets=[0, 1], ancillas1=[2, 3],
                                 ancillas2=[4, 5], sigma="Z")
    double_sel.add_circuit(circuit_block=cb.swap_pair,
                           pair=[0, 1])
    double_sel.add_circuit(circuit_block=epl.append_circuit)
    double_sel.add_circuit(circuit_block=epl.append_circuit)
    double_sel.add_circuit(circuit_block=cb.start_epl)

    return double_sel


def ghz4_epl(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using 4 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)


    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5, 6, 7])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
                    targets=[4, 5, 6, 7], sigma="X")

    # Phase 2 Create last two pairs
    pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                           pair=[0, 1])
    pair.add_circuit(circuit_block=cb.start_epl)
    wrap_EPL_parallel = circuit.Circuit(a0=a0, a1=a1,
                                        circuit_block=pair.run_parallel)

    ghz.add_circuit(circuit_block=wrap_EPL_parallel.append_circuit)

    # Phase 1 Create initial two pairs

    ghz.add_circuit(circuit_block=epl.run_parallel)

    return ghz


def ghz4_epl_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using 3 Bell pairs
    generated using the EPL protocol.
    Note: On simplified you need to substract the average time it
    takes to create a EPL pair

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)


    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
                    targets=[4, 5], sigma="X")

    # Phase 2 Create last two pairs
    pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                           pair=[0, 1])
    pair.add_circuit(circuit_block=cb.start_epl)
    ghz.add_circuit(circuit_block=pair.append_circuit)

    # Phase 1 Create initial two pairs

    ghz.add_circuit(circuit_block=epl.run_parallel)

    return ghz


def ghz4_bk(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using 4 Bell pairs
    generated using the Barret-Kok protocol.

    Parameters
    ----------
    ps : single qubit gate error rate.
    pm : measurement error rate.
    pg : two qubit gate error rate.
    eta : detection efficiency.
    a0 : extra environmental error when electron spin is being operated.
    a1 : default environmental error.
    theta : determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5, 6, 7])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
                    targets=[4, 5, 6, 7], sigma="X")

    # Phase 2 Create last two pairs
    pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                           pair=[0, 1])
    pair.add_circuit(circuit_block=cb.start_BK)
    wrap_EPL_parallel = circuit.Circuit(a0=a0, a1=a1,
                                        circuit_block=pair.run_parallel)

    ghz.add_circuit(circuit_block=wrap_EPL_parallel.append_circuit)

    # Phase 1 Create initial two pairs
    bk = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_BK)
    ghz.add_circuit(circuit_block=bk.run_parallel)

    return ghz


def ghz4_bk_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using 3 Bell pairs
    generated using the Barret-Kok protocol.

    Parameters
    ----------
    ps : single qubit gate error rate.
    pm : measurement error rate.
    pg : two qubit gate error rate.
    eta : detection efficiency.
    a0 : extra environmental error when electron spin is being operated.
    a1 : default environmental error.
    theta : determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
                    targets=[4, 5], sigma="X")

    # Phase 2 Create last two pairs
    pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                           pair=[0, 1])
    pair.add_circuit(circuit_block=cb.start_BK)
    ghz.add_circuit(circuit_block=pair.append_circuit)

    # Phase 1 Create initial two pairs

    bk = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_BK)
    ghz.add_circuit(circuit_block=bk.run_parallel)

    return ghz


def ghz4_single(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using using the single selection
    protocol (see documentation).

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_single_sel(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5, 6, 7])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
                    targets=[4, 5, 6, 7], sigma="X")

    # Phase 2 Create last two pairs
    swaped_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                  pair=[0, 1])
    swaped_pair.add_circuit(circuit_block=pair.run_as_block)
    wrap_pair_parallel = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=swaped_pair.run_parallel)
    ghz.add_circuit(circuit_block=wrap_pair_parallel.append_circuit)

    # Phase 1 Create initial two pairs
    ghz.add_circuit(circuit_block=pair.run_parallel)

    return ghz


def ghz4_single_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using 3 pairs
    using the single selection
    protocol (see documentation).

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_single_sel(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
                    targets=[4, 5], sigma="X")

    # Phase 2 Create last two pairs
    swaped_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                  pair=[0, 1])
    swaped_pair.add_circuit(circuit_block=pair.run_as_block)
    ghz.add_circuit(circuit_block=swaped_pair.append_circuit)

    # Phase 1 Create initial two pairs
    ghz.add_circuit(circuit_block=pair.run_parallel)

    return ghz


def ghz4_double(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using 4 pairs
    using the double selection
    protocol (see documentation).

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_double_sel(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5, 6, 7])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3, 0, 2],
                    targets=[4, 5, 6, 7], sigma="X")

    # Phase 2 Create last two pairs
    swaped_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                  pair=[0, 1])
    swaped_pair.add_circuit(circuit_block=pair.run_as_block)
    wrap_pair_parallel = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=swaped_pair.run_parallel)
    ghz.add_circuit(circuit_block=wrap_pair_parallel.append_circuit)

    # Phase 1 Create initial two pairs
    ghz.add_circuit(circuit_block=pair.run_parallel)

    return ghz


def ghz4_double_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 4 created using 3 pairs
    using the double selection
    protocol (see documentation).

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_double_sel(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=4,
                          measure_pos=[4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 3],
                    targets=[4, 5], sigma="X")

    # Phase 2 Create last two pairs
    swaped_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                  pair=[0, 1])
    swaped_pair.add_circuit(circuit_block=pair.run_as_block)
    ghz.add_circuit(circuit_block=swaped_pair.append_circuit)

    # Phase 1 Create initial two pairs
    ghz.add_circuit(circuit_block=pair.run_parallel)

    return ghz


def ghz3_epl_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 2 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)


    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[2])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
                    targets=[2], sigma="X")

    # Phase 2 Create second pair
    ghz.add_circuit(circuit_block=cb.swap_pair,
                    pair=[0, 1])
    ghz.add_circuit(circuit_block=epl.append_circuit)

    # Phase 1 Create initial pair
    ghz.add_circuit(circuit_block=cb.start_epl)

    return ghz


def ghz3_epl(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 3 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[3, 4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 2, 0],
                    targets=[3, 4, 5], sigma="X")

    # Phase 1 Create initial pair
    start_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                 pair=[0, 1])
    start_pair.add_circuit(circuit_block=cb.start_epl)

    ghz.add_circuit(circuit_block=start_pair.run_parallel3)

    return ghz


def ghz3_bk_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 2 Bell pairs
    generated using the Barret-Kok protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)


    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[2])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
                    targets=[2], sigma="X")

    # Phase 2 Create second pair
    ghz.add_circuit(circuit_block=cb.swap_pair,
                    pair=[0, 1])
    bk = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_BK)
    ghz.add_circuit(circuit_block=bk.append_circuit)

    # Phase 1 Create initial pair
    ghz.add_circuit(circuit_block=cb.start_BK)

    return ghz


def ghz3_bk(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 3 Bell pairs
    generated using the Barret-Kok protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[3, 4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 2, 0],
                    targets=[3, 4, 5], sigma="X")

    # Phase 1 Create initial pair
    start_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                 pair=[0, 1])
    start_pair.add_circuit(circuit_block=cb.start_BK)

    ghz.add_circuit(circuit_block=start_pair.run_parallel3)

    return ghz

def ghz3_single_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 2 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_single_sel(ps, pm, pg, eta, a0, a1, theta)


    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[2])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
                    targets=[2], sigma="X")

    # Phase 2 Create second pair
    ghz.add_circuit(circuit_block=cb.swap_pair,
                    pair=[0, 1])
    ghz.add_circuit(circuit_block=pair.append_circuit)

    # Phase 1 Create initial pair
    ghz.add_circuit(circuit_block=pair.start)

    return ghz


def ghz3_single(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 3 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_single_sel(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[3, 4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 2, 0],
                    targets=[3, 4, 5], sigma="X")

    # Phase 1 Create initial pair
    swaped_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                  pair=[0, 1])
    swaped_pair.add_circuit(circuit_block=pair.run_as_block)
    ghz.add_circuit(circuit_block=swaped_pair.run_parallel3)

    return ghz

def ghz3_double_simple(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 2 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_double_sel(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[2])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
                    targets=[2], sigma="X")

    # Phase 2 Create second pair
    ghz.add_circuit(circuit_block=cb.swap_pair,
                    pair=[0, 1])
    ghz.add_circuit(circuit_block=pair.append_circuit)

    # Phase 1 Create initial pair
    ghz.add_circuit(circuit_block=pair.start)

    return ghz


def ghz3_double(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 3 created using 3 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    # Circuits are assemebled in reversed order
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_double_sel(ps, pm, pg, eta, a0, a1, theta)

    # Phase 3 - Create GHZ
    # Perform the measurements
    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.collapse_ancillas_GHZ,
                          ghz_size=3,
                          measure_pos=[3, 4, 5])
    # Apply two qubit gates in the nodes
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1, 2, 0],
                    targets=[3, 4, 5], sigma="X")

    # Phase 1 Create initial pair
    swaped_pair = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.swap_pair,
                                  pair=[0, 1])
    swaped_pair.add_circuit(circuit_block=pair.run_as_block)
    ghz.add_circuit(circuit_block=swaped_pair.run_parallel3)

    return ghz

def ghz2_epl(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 2 created using 2 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    epl = pair_EPL(ps, pm, pg, eta, a0, a1, theta)

    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.single_selection_ops,
                          targets=[0, 1], ancillas=[2, 3], sigma="Z")
    ghz.add_circuit(circuit_block=cb.swap_pair,
                    pair=[0, 1])
    ghz.add_circuit(circuit_block=epl.run_parallel)

    return ghz

def ghz2_single(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 2 created using 2 Bell pairs
    generated using the EPL protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_single_sel(ps, pm, pg, eta, a0, a1, theta)

    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.single_selection_ops,
                          targets=[0, 1], ancillas=[2, 3], sigma="Z")
    ghz.add_circuit(circuit_block=cb.swap_pair,
                    pair=[0, 1])
    ghz.add_circuit(circuit_block=pair.run_parallel)

    return ghz

def ghz2_double(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 2 created using 2 Bell pairs
    generated using the double selection protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_double_sel(ps, pm, pg, eta, a0, a1, theta)

    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.single_selection_ops,
                          targets=[0, 1], ancillas=[2, 3], sigma="Z")
    ghz.add_circuit(circuit_block=cb.swap_pair,
                    pair=[0, 1])
    ghz.add_circuit(circuit_block=pair.run_parallel)

    return ghz

def ghz2_bk(ps, pm, pg, eta, a0, a1, theta):
    """
    GHZ state of weigth 2 created using 2 Bell pairs
    generated using the double selection protocol.

    Parameters
    ----------
    ps : (scalar) single qubit gate error rate.
    pm : (scalar) measurement error rate.
    pg : (scalar) two qubit gate error rate.
    eta : (scalar) detection efficiency.
    a0 : (scalar) extra environmental error when electron spin is being operated.
    a1 : (scalar) default environmental error.
    theta : (scalar) determines how the states are initialized when generating remote
            entanglement.
    """
    cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
    pair = pair_double_sel(ps, pm, pg, eta, a0, a1, theta)

    ghz = circuit.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_BK)

    return ghz
