"""
Kraus decomposition routines.
"""
import numpy as np
import qutip as qt
import itertools
import operations

error_symbols = ["X", "Y", "Z"]

def char_to_err(c):
    if c == "X":
        return qt.rx
    elif c == "Y":
        return qt.ry
    elif c == "Z":
        return qt.rz
    else:
        return None

def string_to_err(s, positions, system_size):
    s_list = list(s)
    err = qt.qeye([2] * system_size)
    for i in range(len(s_list)):
        e_func = char_to_err(s_list[i])
        if e_func:
            err *= e_func(np.pi, system_size, positions[i])
    return err

def error_combinations(n_errors, err_positions, system_size):
    e_combinations = itertools.combinations_with_replacement(error_symbols, n_errors)
    all_errors = []
    l_errs = len(err_positions)
    for ec in e_combinations:
        err = ''.join(ec)
        permutations = []
        if len(err) < l_errs:
            err = err + "I"*(l_errs - len(err))
        err_permutations = [''.join(p) for p in itertools.permutations(err)]
        err_permutations = set(err_permutations)
        for i in err_permutations:
            single_err = (string_to_err(i, err_positions, system_size), i)
            permutations += [single_err]
        all_errors += [permutations]
    return all_errors


def get_pauli_errors(error_positions, system_size):
    all_errors = [["A", error_combinations(1, error_positions, system_size)],
                  ["B", error_combinations(2, error_positions, system_size)],
                  ["C", error_combinations(3, error_positions, system_size)],
                  ["D", error_combinations(4, error_positions, system_size)]]
    return all_errors
