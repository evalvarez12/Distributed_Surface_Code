"""
Kraus decomposition routines.
"""
import numpy as np
import qutip as qt
import operations

sigmas = [qt.rx, qt.ry, qt.rz]
symbols = ["X", "Y", "Z"]
L = len(sigmas)
# def operational_state(N):
#     state = qt.tensor(qt.basis(2**N, 0), qt.basis(2**N, 0))
#     for i in range(1, 2**N):
#         state += qt.tensor(qt.basis(2**N, i), qt.basis(2**N, i))
#     return 1/np.sqrt(2**N) * state


def get_pauli_errors(system_size, acting_positions):
    all_errors = [error_A(system_size, acting_positions), error_B(system_size, acting_positions),
                  error_C(system_size, acting_positions), error_D(system_size, acting_positions)]
    return all_errors

def error_A(system_size, acting_positions):
    """
    Error A corresponds to single qubit error
    """
    L_pos = len(acting_positions)

    # List of all possible errors
    all_errors = []
    error_names = []
    # Choose a error
    for err_ind in range(L):
        # Find all possible permutations
        error_permutations = []
        for pos_ind in range(L_pos):
            error_permutations += [sigmas[err_ind](np.pi, system_size, acting_positions[pos_ind])]
        all_errors += [error_permutations]
        error_names += [symbols[err_ind]]

    return all_errors, error_names

def error_B(system_size, acting_positions):
    """
    Error B corresponds to two qubits error
    """
    L_pos = len(acting_positions)

    # List of all possible errors
    all_errors = []
    error_names = []
    # Choose a error
    for err1_ind in range(L):
        for err2_ind in range(L):
            # Find all possible permutations
            error_permutations = []
            for pos1_ind in range(L_pos):
                for pos2_ind in range(pos1_ind + 1, L_pos):
                    error_permutations += [sigmas[err1_ind](np.pi, system_size, acting_positions[pos1_ind])
                                           * sigmas[err2_ind](np.pi, system_size, acting_positions[pos2_ind])]
            all_errors += [error_permutations]
            error_names += [symbols[err1_ind] + symbols[err2_ind]]

    return all_errors, error_names


def error_C(system_size, acting_positions):
    """
    Error C corresponds to three qubits error
    """
    L_pos = len(acting_positions)

    # List of all possible errors
    all_errors = []
    error_names = []
    # Choose a error
    for err1_ind in range(L):
        for err2_ind in range(L):
            for err3_ind in range(L):
                # Find all possible permutations
                error_permutations = []
                for pos1_ind in range(L_pos):
                    for pos2_ind in range(pos1_ind + 1, L_pos):
                        for pos3_ind in range(pos2_ind + 1, L_pos):
                            error_permutations += [sigmas[err1_ind](np.pi, system_size, acting_positions[pos1_ind])
                                                   * sigmas[err2_ind](np.pi, system_size, acting_positions[pos2_ind])
                                                   * sigmas[err3_ind](np.pi, system_size, acting_positions[pos3_ind])]
                all_errors += [error_permutations]
                error_names += [symbols[err1_ind] + symbols[err2_ind]
                              + symbols[err3_ind]]

    return all_errors, error_names


def error_D(system_size, acting_positions):
    """
    Error D corresponds to four qubits error
    """
    L_pos = len(acting_positions)

    # List of all possible errors
    all_errors = []
    error_names = []
    # Choose a error
    for err1_ind in range(L):
        for err2_ind in range(L):
            for err3_ind in range(L):
                for err4_ind in range(L):
                    # Find all possible permutations
                    error_permutations = []
                    for pos1_ind in range(L_pos):
                        for pos2_ind in range(pos1_ind + 1, L_pos):
                            for pos3_ind in range(pos2_ind + 1, L_pos):
                                for pos4_ind in range(pos3_ind + 1, L_pos):
                                    error_permutations += [sigmas[err1_ind](np.pi, system_size, acting_positions[pos1_ind])
                                                           * sigmas[err2_ind](np.pi, system_size, acting_positions[pos2_ind])
                                                           * sigmas[err3_ind](np.pi, system_size, acting_positions[pos3_ind])
                                                           * sigmas[err4_ind](np.pi, system_size, acting_positions[pos4_ind])]
                    all_errors += [error_permutations]
                    error_names += [symbols[err1_ind] + symbols[err2_ind]
                                  + symbols[err4_ind] + symbols[err4_ind]]

    return all_errors, error_names
