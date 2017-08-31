"""
Routines to perform entanglement purification.

author: Eduardo Villaseñor
created-on: 05/06/17
"""

import qutip as qt
import numpy as np


def random_measure_single_Xbasis(rho, N=1, pos=0, dimRed=False):
    """
    Measure a single qubit in the X basis.

    Parameters
    ----------
    rho : density matrix.
    N : system size.
    pos : position of qubit to be measured.
    dimRed : is the collapsed state has reduced to the
            dimentions N - 1
    """
    H = qt.snot(N, pos)
    rho = H * rho * H.dag()
    measurement, collapsed_rho = random_measure_single_Zbasis(rho, N, pos, dimRed)
    if not dimRed:
        collapsed_rho = H * rho * H.dag()
    return measurement, collapsed_rho

def random_measure_single_Zbasis(rho, N=1, pos=0, dimRed=False):
    """
    Measure a single qubit in the Z basis, using a random number to
    choose a measurement outcome.

    Parameters
    ----------
    rho : density matrix.
    N : system size.
    pos : position of qubit to be measured.
    dimRed : is the collapsed state has reduced to the
            dimentions N - 1
    """
    # Calculate the probability of measurment 1
    p0 = p_measurement_single_Zbasis(rho, 0, N, pos, dimRed)
    r = np.random.rand()
    # print("P0", p0)
    # Draw a measurement
    if r < np.linalg.norm(p0):
        collapsed_rho = collapse_single_Zbasis(rho, 0, N, pos, dimRed)
        measurement = 1
    else:
        collapsed_rho = collapse_single_Zbasis(rho, 1, N, pos, dimRed)
        measurement = -1
    return measurement, collapsed_rho


def collapse_single_Zbasis(rho, proyect, N=1, pos=0, dimRed=False):
    """
    Collapse the state in the postion of a single qubit.
    """
    # Obtain the proyection operator depending on dimRed
    if dimRed:
        p = proyector_single_qubit_Zbasis_dimRed(proyect, N, pos)
    else:
        p = proyector_single_qubit_Zbasis(proyect, N, pos)
    collapsed_rho = p * rho * p.dag()
    return collapsed_rho/collapsed_rho.tr()


def p_measurement_single_Zbasis(rho, measure, N=1, pos=0, dimRed=False):
    """
    Calculate the probability of measuring the value "measure".
    """
    p = proyector_single_qubit_Zbasis(measure, N, pos)
    return (p * rho).tr()


def proyector_single_qubit_Zbasis_dimRed(proyect, N=1, pos=0):
    # NOTE: This proyector reduces the dimension of the state density matrix
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
    """
    Tensor a single qubit operator between identyties accoring to the
    position of the qubit and the system size.
    """
    # TODO use qt.rx, qt.ry, qt.rz instead for sigmas?
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
    """
    Tensor the product of multiple single qubit operators.
    """
    l = len(operators)
    res = qt.qeye([2]*N)
    for i in range(l):
        res *= tensor_single_operator(operators[i], N, positions[i])
    return res
