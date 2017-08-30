"""
Kraus decomposition routines.
"""
import numpy as np
import qutip as qt
import operations

sigmas = [qt.rx, qt.ry, qt.rz]


# def operational_state(N):
#     state = qt.tensor(qt.basis(2**N, 0), qt.basis(2**N, 0))
#     for i in range(1, 2**N):
#         state += qt.tensor(qt.basis(2**N, i), qt.basis(2**N, i))
#     return 1/np.sqrt(2**N) * state


def error_A(system_size):
    """
    Error A corresponds to single qubit error
    """
    # List of all possible errors
    all_errors = []
    # Choose a error
    for err in sigmas:
        # Find all possible permutations
        error_permutations = []
        for pos in range(system_size):
            error_permutations += [err(np.pi, system_size, pos)]
        all_errors += [error_permutations]

    return all_errors

def error_B(system_size):
    """
    Error B corresponds to two qubits error
    """
    # List of all possible errors
    all_errors = []
    # Choose a error
    for err1 in sigmas:
        for err2 in sigmas:
            # Find all possible permutations
            error_permutations = []
            for pos1 in range(system_size):
                for pos2 in range(pos1 + 1, system_size):
                    error_permutations += [err1(np.pi, system_size, pos1)
                                           * err2(np.pi, system_size, pos2)]
            all_errors += [error_permutations]

    return all_errors


def error_C(system_size):
    """
    Error C corresponds to three qubits error
    """
    # List of all possible errors
    all_errors = []
    # Choose a error
    for err1 in sigmas:
        for err2 in sigmas:
            for err3 in sigmas:
                # Find all possible permutations
                error_permutations = []
                for pos1 in range(system_size):
                    for pos2 in range(pos1 + 1, system_size):
                        for pos3 in range(pos2 + 1, system_size):
                            error_permutations += [err1(np.pi, system_size, pos1)
                                                   * err2(np.pi, system_size, pos2)
                                                   * err3(np.pi, system_size, pos3)]
                all_errors += [error_permutations]

    return all_errors


def error_D(system_size):
    """
    Error D corresponds to four qubits error
    """
    # List of all possible errors
    all_errors = []
    # Choose a error
    for err1 in sigmas:
        for err2 in sigmas:
            for err3 in sigmas:
                for err4 in sigmas:
                    # Find all possible permutations
                    error_permutations = []
                    for pos1 in range(system_size):
                        for pos2 in range(pos1 + 1, system_size):
                            for pos3 in range(pos2 + 1, system_size):
                                for pos4 in range(pos3 + 1, system_size):
                                    error_permutations += [err1(np.pi, system_size, pos1)
                                                           * err2(np.pi, system_size, pos2)
                                                           * err3(np.pi, system_size, pos3)
                                                           * err4(np.pi, system_size, pos4)]
                    all_errors += [error_permutations]

    return all_errors


def fidelity(rhoA, stateB):
    """
    Fidelity for the special case when one of the states is a pure state.
    """
    return stateB.dag() * rhoA * stateB
