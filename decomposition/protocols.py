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
import protocols_det

# Determine parameters
ps = 0.003
pm = 0.003
pg = 0.003
# a0 = 83.33
a0 = 8.0
# a1 = 1/3.
a1 = 1/80.
# eta = (0.1)*(0.03)*(0.8)
eta = 1/100.
theta = .86

iterations = 50

# Initialize objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
# cb = circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4.)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
rho_ref2 = qt.bell_state('01') * qt.bell_state('01').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()
ghz3_ref = qt.ghz_state(3) * qt.ghz_state(3).dag()
pd = protocols_det.ProtocolsDeterministic(0,0,0,0)

# Lists to save results
# SIMPLE_fidelity = []
# SIMPLE_steps = []
# MEDIUM_fidelity = []
# MEDIUM_steps = []
# COMPLEX_fidelity = []
# COMPLEX_steps = []
FIDELITY = []
TIMES = []

# for eta in [1/30., 1/40., 1/50., 1/60., 1/70., 1/80.]:
# for a0 in [12., 10., 8., 6., 4., 2.]:
# for a0 in [40., 30., 20., 10., 5., 2.]:
for nada in [1]:
    cb.change_parameters(ps, pm, pg, eta, a0, a1, theta)
    print("-------------------var =" + str(eta) + "------------------")

    '''PROTOCOLS START HERE'''
    # print("-------------------PROTOCOL SIMPLEST BK 3------------------")
    # # First assemeble the small single selection circuit
    # start_BK = circuit.Circuit(a0=a0, a1=a1,
    #                                      circuit_block=cb.start_BK)
    # start_BK.add_circuit(circuit_block=cb.swap_pair,
    #                                pair=[0, 1])
    #
    #
    # # Phase 1 - Purify Bell pair
    # ghz = circuit.Circuit(a0=a0, a1=a1,
    #                       circuit_block=cb.start_BK)
    #
    # # Phase 2 - Create GHZ
    # ghz.add_circuit(circuit_block=start_BK.append_circuit)
    # ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[1],
    #                 targets=[2], sigma="X")
    # ghz.add_circuit(circuit_block=cb.collapse_ancillas_Z,
    #                 ancillas_pos=[2],
    #                 projections=[0])
    #
    # # Get average number of steps
    # fidelity = []
    # times = []
    # rho = ghz3_ref*0
    # # check = collections.Counter({})
    # for i in range(iterations):
    #     _, c, r = ghz.run(None)
    #     times += [c["time"]]
    #
    #     r = pd.twirl_ghz(r)
    #     fidelity += [qt.fidelity(r, ghz3_ref)]
    #
    #     rho += r
    #
    # print("F: ", np.average(fidelity), np.std(fidelity))
    # print("TIMES:", np.average(times))
    # FIDELITY += [(np.average(fidelity), np.std(fidelity))]
    # TIMES += [(np.average(times), np.std(times))]
    # name = "ghz_3_eta_" + str(round(eta, 3))
    # # qt.qsave(rho, name)

# print(rho)
# SIMPLE_fidelity += [(np.average(fidelity), np.std(fidelity))]
# SIMPLE_steps += [(np.average(steps), np.std(steps))]
# print("End protocol")
# FIDELITY = np.array(FIDELITY)
# TIMES = np.array(TIMES)
#
# np.save("3eta_thres_fidelity.npy", FIDELITY)
# np.save("3eta_thres_times.npy", TIMES)
# Average over values
#



    '''PROTOCOLS OF GHZ 4'''
    # print("-------------------PROTOCOL SIMPLEST BK------------------")
    # # First assemeble the small single selection circuit
    # single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
    #                                      circuit_block=cb.start_BK)
    # single_sel_simple2.add_circuit(circuit_block=cb.swap_pair,
    #                                pair=[0, 1])
    # wrap_single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
    #                                           circuit_block=single_sel_simple2.run_parallel)
    #
    # single_sel_simple1 = circuit.Circuit(a0=a0, a1=a1,
    #                                      circuit_block=cb.start_BK)
    #
    # # Phase 1 - Purify Bell pair
    # ghz = circuit.Circuit(a0=a0, a1=a1,
    #                       circuit_block=single_sel_simple1.run_parallel)
    #
    # # Phase 2 - Create GHZ
    # ghz.add_circuit(circuit_block=wrap_single_sel_simple2.append_circuit)
    # ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
    #                 targets=[1, 3, 0, 2], sigma="Z")
    # ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
    #                 ancillas_pos=[4, 5, 6, 7],
    #                 projections=[0, 0, 0, 0])
    #
    # # Get average number of steps
    # fidelity = []
    # times = []
    # rho = ghz_ref*0
    # # check = collections.Counter({})
    # for i in range(iterations):
    #     _, c, r = ghz.run(None)
    #     times += [c["time"]]
    #
    #     r = pd.twirl_ghz(r)
    #     fidelity += [qt.fidelity(r, ghz_ref)]
    #
    #     rho += r
    #
    # print("F: ", np.average(fidelity), np.std(fidelity))
    # print("TIMES:", np.average(times))
    # FIDELITY += [(np.average(fidelity), np.std(fidelity))]
    # TIMES += [(np.average(times), np.std(times))]
    # name = "ghz_4_a0_" + str(round(a0, 3))
    # # qt.qsave(rho, name)

# print(rho)
# SIMPLE_fidelity += [(np.average(fidelity), np.std(fidelity))]
# SIMPLE_steps += [(np.average(steps), np.std(steps))]
# print("End protocol")
# FIDELITY = np.array(FIDELITY)
# TIMES = np.array(TIMES)
#
# np.save("a_thres_fidelity.npy", FIDELITY)
# np.save("a_thres_times.npy", TIMES)
# Average over values

    print("-------------------PROTOCOL SIMPLEST EPL------------------")
    # First assemeble the small single selection circuit
    single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=cb.start_epl)

    wrap_single_sel_simple2 = circuit.Circuit(a0=a0, a1=a1,
                                              circuit_block=single_sel_simple2.run_parallel)

    single_sel_simple1 = circuit.Circuit(a0=a0, a1=a1,
                                         circuit_block=cb.start_epl)


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

    # Get average number of steps
    fidelity = []
    check = collections.Counter({})
    for i in range(iterations):
        p, c, rho = ghz.run(None)
        check += c
        fidelity += [qt.fidelity(rho, ghz_ref)]

    # print(rho)
    # SIMPLE_fidelity += [(np.average(fidelity), np.std(fidelity))]
    # SIMPLE_steps += [(np.average(steps), np.std(steps))]
    print("End protocol")

    # Average over check values
    for k in check:
        check[k] = check[k]/iterations
    print("TIME: ", check["time"])
    print("F: ", np.average(fidelity), np.std(fidelity))


    print("-------------------PROTOCOL SIMPLE------------------")
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

    # Get average number of steps
    fidelity = []
    check = collections.Counter({})
    for i in range(iterations):
        p, c, rho = ghz.run(None)
        check += c
        fidelity += [qt.fidelity(rho, ghz_ref)]

    # print(rho)
    # SIMPLE_fidelity += [(np.average(fidelity), np.std(fidelity))]
    # SIMPLE_steps += [(np.average(steps), np.std(steps))]
    print("End protocol")

    # Average over check values
    for k in check:
        check[k] = check[k]/iterations
    print("TIME: ", check["time"])
    print("F: ", np.average(fidelity), np.std(fidelity))


    """
    Nickerson expedient protocol.
    """
    print("------------------PROTOCOL MEDIUM-------------------")
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
    # MEDIUM_fidelity += [(np.average(fidelity), np.std(fidelity))]
    # MEDIUM_steps += [(np.average(steps), np.std(steps))]
    print("End protocol")


    print("TIME: ", check["time"])
    print("F: ", np.average(fidelity), np.std(fidelity))


    """
    Nickerson stringent protocol.
    """
    print("------------------PROTOCOL COMPLEX-------------------")
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

    # Phase 1 - Purify Bell pair
    ghz = circuit.Circuit(a0=a0, a1=a1, circuit_block=double_sel1.run_parallel)


    # Phase 2 - Create GHZ
    ghz.add_circuit(circuit_block=wrap_single_sel_medium2.append_circuit)
    ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
                    targets=[1, 3, 0, 2], sigma="Z")

    ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
                    ancillas_pos=[4, 5, 6, 7],
                    projections=[0, 0, 0, 0])
    # ghz.add_circuit(circuit_block=wrap_single_sel_medium2.append_circuit)
    # ghz.add_circuit(circuit_block=cb.two_qubit_gates, controls=[4, 5, 6, 7],
    #                 targets=[1, 3, 0, 2], sigma="Z")
    #
    # ghz.add_circuit(circuit_block=cb.collapse_ancillas_X,
    #                 ancillas_pos=[4, 5, 6, 7],
    #                 projections=[0, 0, 0, 0])

    # Get average number of steps
    fidelity = []
    check = collections.Counter({})
    for i in range(iterations):
        p, c, rho = ghz.run(None)
        check += c
        fidelity += [qt.fidelity(rho, ghz_ref)]

    # COMPLEX_fidelity += [(np.average(fidelity), np.std(fidelity))]

    for k in check:
        check[k] = check[k]/iterations

    # COMPLEX_steps += [(np.average(steps), np.std(steps))]
    # print(rho)
    print("End protocol")
    print("TIME: ", check["time"])
    print("F: ", np.average(fidelity), np.std(fidelity))

# np.save("data/SIMPLE_fidelity_BK", SIMPLE_fidelity)
# np.save("data/SIMPLE_steps_BK", SIMPLE_steps)
# np.save("data/MEDIUM_fidelity_BK", MEDIUM_fidelity)
# np.save("data/MEDIUM_steps_BK", MEDIUM_steps)
# np.save("data/COMPLEX_fidelity_BK", COMPLEX_fidelity)
# np.save("data/COMPLEX_steps_BK", COMPLEX_steps)
# np.save("data/p_env", p_env_var)
