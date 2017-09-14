"""
Kraus decomposition routines.
"""
import numpy as np
import qutip as qt
import itertools
import operations

error_symbols = ["X", "Y", "Z"]

def char_to_err(c):
    """Converts a character into a sigma function."""
    if c == "X":
        return qt.rx
    elif c == "Y":
        return qt.ry
    elif c == "Z":
        return qt.rz
    else:
        return None

def string_to_err(s, positions, system_size):
    """
    Takes a string and returns a gate acting on the corresponding qubits.
    Strings example: IXZI, YYYZ, ...
    """
    s_list = list(s)
    err = qt.qeye([2] * system_size)
    for i in range(len(s_list)):
        e_func = char_to_err(s_list[i])
        if e_func:
            err *= e_func(np.pi, system_size, positions[i])
    return err

def error_combinations(n_errors, err_positions, system_size):
    """
    Returns all the possible error combinations for a given number of errors.
    """
    l_errs = len(err_positions)

    # If no errors generate the identity manually
    if n_errors == 0:
        return [[(qt.qeye([2] * system_size), "I" * l_errs)]]

    # Get all the possible combinations of errors without permutations
    e_combinations = itertools.combinations_with_replacement(error_symbols,
                                                             n_errors)
    # List to save all errors
    all_errors = []
    for ec in e_combinations:
        # Join the selected error into a single string
        err = ''.join(ec)
        # Sub list to save different permutations
        permutations = []
        # Append the extra needed identity symbols
        if len(err) < l_errs:
            err = err + "I"*(l_errs - len(err))
        # Find the permutations of a selected error and remove the duplicates
        err_permutations = [''.join(p) for p in itertools.permutations(err)]
        err_permutations = set(err_permutations)
        # Get the error gate for each case
        for i in err_permutations:
            # Both error and symbol are saved
            single_err = (string_to_err(i, err_positions, system_size), i)
            permutations += [single_err]
        all_errors += [permutations]
    return all_errors


def get_pauli_errors(error_positions, system_size):
    """Get all errors from 0 to 4."""
    all_errors = [[error_combinations(0, error_positions, system_size), "I"],
                  [error_combinations(1, error_positions, system_size), "A"],
                  [error_combinations(2, error_positions, system_size), "B"],
                  [error_combinations(3, error_positions, system_size), "C"],
                  [error_combinations(4, error_positions, system_size), "D"]]
    return all_errors

def symbol_product(symbolA, symbolB):
    # Turn symbols into lists
    list_A = list(symbolA)
    list_B = list(symbolB)
    prod = list_A.copy()
    l = len(list_A)

    # Apply product for each element in list
    for i in range(l):
        prod[i] = product_rules(list_A[i], list_B[i])

    # Return the result in a single string
    return ''.join(prod)

def product_rules(A, B):
    # Check first for identities
    if A == "I":
        return B
    if B == "I":
        return A

    # Check arguments
    if not A in error_symbols or not B in error_symbols:
        raise ValueError("Arguments are not in {X, Y, Z}")

    # Check if they are equal
    if A == B:
        return "I"

    # Find the missing element and return it
    prod = [A, B]
    res = list(set(error_symbols) - set(prod))[0]
    return res
