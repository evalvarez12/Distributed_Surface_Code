"""
Routines to assemble the pauli basis based on constructiong the
operators based on symbols ("IIXI", "YZZX", ... )
"""
import numpy as np
import qutip as qt
import itertools
import operations

operator_symbols = ["X", "Y", "Z"]

def char_to_operator(c):
    """Converts a character into a sigma function."""
    if c == "X":
        return qt.rx
    elif c == "Y":
        return qt.ry
    elif c == "Z":
        return qt.rz
    else:
        return None

def string_to_operator(s, positions, system_size):
    """
    Takes a string and returns a gate acting on the corresponding qubits.
    Strings example: IXZI, YYYZ, ...
    """
    s_list = list(s)
    operator = qt.qeye([2] * system_size)
    for i in range(len(s_list)):
        op_func = char_to_operator(s_list[i])
        if op_func:
            operator *= op_func(np.pi, system_size, positions[i])
    return operator

def combinations(n_operators, ops_positions, system_size):
    """
    Returns all the possible operator combinations for a given number of operators.
    """
    # Dict to save all operators
    operators = {}
    l_ops = len(ops_positions)

    # If no operators generate the identity manually
    if n_operators == 0:
        operators["I" * l_ops] = qt.qeye([2] * system_size)
        return operators

    # Get all the possible combinations of operators without permutations
    e_combinations = itertools.combinations_with_replacement(operator_symbols,
                                                             n_operators)
    for ec in e_combinations:
        # Join the selected operator into a single string
        ops = ''.join(ec)

        # Append the extra needed identity symbols
        if len(ops) < l_ops:
            ops = ops + "I"*(l_ops - len(ops))
        # Find the permutations of a selected operator and remove the duplicates
        ops_permutations = [''.join(p) for p in itertools.permutations(ops)]
        ops_permutations = set(ops_permutations)
        # Get the operator gate for each case
        for i in ops_permutations:
            # Both operator and symbol are saved
            operators[i] = string_to_operator(i, ops_positions, system_size)
    return operators


def get_basis(operator_positions, system_size):
    """Get all operators from 0 to system_size."""
    # Initialize dictionary with the identity
    all_operators = combinations(0, operator_positions, system_size)
    # Update dict with all combinations for each number of identities
    for i in range(1, len(operator_positions) + 1):
        all_operators.update(combinations(i, operator_positions, system_size))
    return all_operators

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
    if not A in operator_symbols or not B in operator_symbols:
        raise ValueError("Arguments are not in {X, Y, Z}")

    # Check if they are equal
    if A == B:
        return "I"

    # Find the missing element and return it
    prod = [A, B]
    res = list(set(operator_symbols) - set(prod))[0]
    return res
