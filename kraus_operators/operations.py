"""
Routines to perform entanglement purification.

author: Eduardo Villaseñor
created-on: 05/06/17
"""

import qutip as qt
import numpy as np


def measure_single_Zbasis(state, N=1, pos=0, dimRed=False):
    p0 = p_measurement_single_Zbasis(state, 0, N, pos, dimRed)
    r = np.random.rand()
    if r < np.linalg.norm(p0):
        measured_state = collapse_single_Zbasis(state, 0, N, pos, dimRed)
        measurement = 1
    else:
        measured_state = collapse_single_Zbasis(state, 1, N, pos, dimRed)
        measurement = -1
    return measurement, measured_state


def collapse_single_Zbasis(state, proyect, N=1, pos=0, dimRed=False):
    if dimRed:
        p = proyector_single_qubit_Zbasis_dimRed(proyect, N, pos)
    else:
        p = proyector_single_qubit_Zbasis(proyect, N, pos)
    collpsed_state = p * state * p.dag()
    return collpsed_state/collpsed_state.tr()


def p_measurement_single_Zbasis(state, measure, N=1, pos=0, dimRed=False):
    # TODO check that if state is not type operator to raise err
    if dimRed:
        p = proyector_single_qubit_Zbasis_dimRed(measure, N, pos)
    else:
        p = proyector_single_qubit_Zbasis(measure, N, pos)
    return (p * state).tr()


def proyector_single_qubit_Zbasis_dimRed(proyect, N=1, pos=0):
    # NOTE: This proyector reduces the dimension of the state
    # is a rectangular matrix <x| not |x><x|
    if proyect != 0 and proyect != 1:
        raise ValueError("proyector: measurement value invalid")

    p = qt.basis(2, proyect).dag()
    return tensor_single_operator(p, N, pos)


def proyector_single_qubit_Zbasis(proyect, N=1, pos=0):
    # NOTE: This proyector is |x><x|
    if proyect != 0 and proyect != 1:
        raise ValueError("proyector: measurement value invalid")

    p = qt.basis(2, proyect)
    p = p * p.dag()
    return tensor_single_operator(p, N, pos)


def tensor_single_operator(operator, N, pos):
    if pos > N:
        raise ValueError("tensor_single_operator: N > pos")
    if N == 1:
        return operator

    if pos == 0:
        idd = qt.qeye([2]*(N-1))
        return qt.tensor(operator, idd)
    elif pos == N - 1:
        idd = qt.qeye([2]*(N-1))
        return qt.tensor(idd, operator)
    else:
        idd1 = qt.qeye([2]*pos)
        idd2 = qt.qeye([2]*(N - pos - 1))
        temp = qt.tensor(idd1, operator)
        return qt.tensor(temp, idd2)


def tensor_operator(operators, positions, N):
    l = len(operators)
    res = qt.qeye([2]*N)
    for i in range(l):
        res *= tensor_single_operator(operators[i], N, positions[i])
    return res
