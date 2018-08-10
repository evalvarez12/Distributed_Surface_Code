"""
Routines to assemble the pauli basis based on constructiong the
operators based on symbols ("IIXI", "YZZX", ... )

author: Eduardo Villaseñor
created-on: 18/11/17
"""
import numpy as np
import qutip as qt
import itertools

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

    Parameters
    ----------
    s : (string) string representation of single qubit operators.
    positions : (list) position in which each operator acts.
    system_size : (int) total system size.
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
    Returns all the possible operator combinations for a given number
    of operators and system size

    Parameters
    ----------
    n_operators : (int) number of operators in which the combinations are going to be
                  obtained.
    ops_positions : (list) position in which each operator acts.
    system_size : (int) total system size.
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
            if i not in operators:
                operators[i] = string_to_operator(i, ops_positions, system_size)
    return operators


def get_basis(operator_positions, system_size):
    """Get all operators comninations from 0 to system_size."""
    # Initialize dictionary with the identity
    all_operators = combinations(0, operator_positions, system_size)
    # Update dict with all combinations for each number of identities
    for i in range(1, len(operator_positions) + 1):
        all_operators.update(combinations(i, operator_positions, system_size))
    return all_operators

def symbol_product(symbolA, symbolB):
    """
    Do the operation of two operators in symbol form up to a -1 factor.
    ex. XYZZ * IXZY = XZIX
    """
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
    """Rulles that used in the symbol product function."""
    # Check first for identities
    if A == "I":
        return B
    if B == "I":
        return A

    # Check arguments
    if A not in operator_symbols or B not in operator_symbols:
        raise ValueError("Arguments are not in {X, Y, Z}")

    # Check if they are equal
    if A == B:
        return "I"

    # Find the missing element and return it
    prod = [A, B]
    res = list(set(operator_symbols) - set(prod))[0]
    return res
