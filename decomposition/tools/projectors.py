"""
Some tools to generate projectors based on even/odd or binary numbers.

author: Eduardo Villaseñor
created-on: 18/11/17
"""
import qutip as qt


def number_to_ket(system_size, number):
    # Returns a qubit ket based on the binary representation of a number
    # Obtain int list binary representation of the number
    bin_list = [int(x) for x in bin(number)[2:]]
    # Attach missing zeros
    if len(bin_list) < system_size:
        bin_list = [0]*(system_size - len(bin_list)) + bin_list

    # NOTE: qutip has the computational basis reversed.
    # ej: rx(pi, 2, 0)  = X tensor Id
    # Create state
    state = qt.basis(2, bin_list[0])
    for i in bin_list[1:]:
        # State is build last qubit to first qubit
        state = qt.tensor(state, qt.basis(2, i))
    return state

def numbers_to_projectors(system_size, numbers):
    # Takes a list of number and returns a list with the projectors of each one|
    ps = []
    for i in numbers:
        ket = number_to_ket(system_size, i)
        ps += [ket * ket.dag()]
    return ps

def numbers_to_kets(system_size, numbers):
    # Returns a list of the ket representation of each number
    kets = []
    for i in numbers:
        ket = number_to_ket(system_size, i)
        kets += [ket]
    return kets

def odd_projectors(system_size):
    """
    Returns all the odd projectors in the qubit basis of a given system size.
    """
    # All the odd numbers in Hamming weight
    odds = []
    for i in range(2**system_size):
        bin_list = [int(x) for x in bin(i)[2:]]
        if bin_list.count(1) % 2 == 1:
            odds += [i]
    odd_projectors = numbers_to_projectors(system_size, odds)
    return sum(odd_projectors)


def even_projectors(system_size):
    """
    Returns all the even projectors in the qubit basis of a given system size.
    """
    # All the even numbers in Hamming weigth
    evens = []
    for i in range(2**system_size):
        bin_list = [int(x) for x in bin(i)[2:]]
        if bin_list.count(1) % 2 == 0:
            evens += [i]
    even_projectors = numbers_to_projectors(system_size, evens)
    return sum(even_projectors)


def project_odd(rho):
    # Project state rho into odd subspace
    system_size = len(rho.dims[0])
    projectors = odd_projectors(system_size)
    projected_rho = rho * 0
    for i in projectors:
        for j in projectors:
            projected_rho += i * rho * j.dag()
    if projected_rho.tr() != 0:
        projected_rho = projected_rho / projected_rho.tr()
    return projected_rho

def project_even(rho):
    # Project state rho into even subspace
    system_size = len(rho.dims[0])
    projectors = even_projectors(system_size)
    projected_rho = rho * 0
    for i in projectors:
        for j in projectors:
            #  P * rho * P.dag()
            projected_rho += i * rho * j
    if projected_rho.tr() != 0:
        projected_rho = projected_rho / projected_rho.tr()
    return projected_rho

def project_odd_ket(state):
    system_size = len(state.dims[0])
    projectors = odd_projectors(system_size)
    projected_state = state * 0
    for i in projectors:
        projected_state += i * state
    if projected_state.norm() != 0:
        projected_state = projected_state / projected_state.norm()
    return projected_state

def project_even_ket(state):
    system_size = len(state.dims[0])
    projectors = even_projectors(system_size)
    projected_state = state * 0
    for i in projectors:
        projected_state += i * state
    if projected_state.norm() != 0:
        projected_state = projected_state / projected_state.norm()
    return projected_state
