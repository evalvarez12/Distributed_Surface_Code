import qutip as qt
import numpy as np


def number_to_ket(system_size, number):
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


def list_numbers_to_projectors(system_size, numbers):
    ps = []
    for i in numbers:
        ket = number_to_ket(system_size, i)
        ps += [ket * ket.dag()]
    return ps


def odd_projectors(system_size):
    # Odd numbers in binary wise
    odds = []
    for i in range(2**system_size):
        bin_list = [int(x) for x in bin(i)[2:]]
        if bin_list.count(1) % 2 == 1:
            odds += [i]
    projectors = list_numbers_to_projectors(system_size, odds)
    return projectors


def even_projectors(system_size):
    # Even numbers in binary wise
    evens = []
    for i in range(2**system_size):
        bin_list = [int(x) for x in bin(i)[2:]]
        if bin_list.count(1) % 2 == 0:
            evens += [i]

    projectors = list_numbers_to_projectors(system_size, evens)
    return projectors


def project_odd(rho):
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
