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


def forced_measure_single_Xbasis(rho, N=1, pos=0, projector=0, dimRed=False):
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
    p, collapsed_rho = forced_measure_single_Zbasis(rho, N, pos, projector, dimRed)
    if not dimRed:
        collapsed_rho = H * rho * H.dag()
    return p, collapsed_rho


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
    # print("P0", p0, r)
    # Draw a measurement
    if r < p0:
        collapsed_rho = collapse_single_Zbasis(rho, 0, N, pos, dimRed)
        measurement = 1
    else:
        collapsed_rho = collapse_single_Zbasis(rho, 1, N, pos, dimRed)
        measurement = -1
    return measurement, collapsed_rho


def forced_measure_single_Zbasis(rho, N=1, pos=0, project=0, dimRed=False):
    """
    Force a measurement on a single qubit in the Z basis to a given proyetor.

    Parameters
    ----------
    rho : density matrix.
    N : system size.
    pos : position of qubit to be measured.
    project : force the measurement into the desired projector
    dimRed : is the collapsed state has reduced to the
            dimentions N - 1
    """
    # Calculate the actual probability of doing the measurement
    p = p_measurement_single_Zbasis(rho, 0, N, pos, dimRed)
    # print("P", project, p)
    if p == 0:
        raise ZeroDivisionError("p = 0!")
    # Draw the measurement
    collapsed_rho = collapse_single_Zbasis(rho, project, N, pos, dimRed)
    return p, collapsed_rho



def collapse_single_Zbasis_ket(psi, project, N=1, pos=0, dimRed=False):
    """
    Collapse the state in the postion of a single qubit in Z basis.
    """
    # Obtain the projection operator depending on dimRed
    if dimRed:
        p = projector_single_qubit_Zbasis_dimRed(project, N, pos)
    else:
        p = projector_single_qubit_Zbasis(project, N, pos)
    collapsed_psi = p * psi
    if collapsed_psi.norm() == 0:
        return collapsed_psi
    return collapsed_psi/collapsed_psi.norm()

def collapse_single_Xbasis_ket(psi, project, N=1, pos=0, dimRed=False):
    """
    Collapse the state in the postion of a single qubit in X basis.
    """
    H = qt.snot(N, pos)
    collapsed_psi = H * psi
    collapsed_psi = collapse_single_Zbasis_ket(collapsed_psi, project, N, pos,
                                               dimRed)
    if not dimRed:
        collapsed_psi = H * collapsed_psi
    return collapsed_psi

def collapse_single_Zbasis(rho, project, N=1, pos=0, dimRed=False):
    """
    Collapse the state in the postion of a single qubit.
    """
    # Obtain the projection operator depending on dimRed
    if dimRed:
        p = projector_single_qubit_Zbasis_dimRed(project, N, pos)
    else:
        p = projector_single_qubit_Zbasis(project, N, pos)
    collapsed_rho = p * rho * p.dag()
    return collapsed_rho/collapsed_rho.tr()


def p_measurement_single_Zbasis(rho, measure, N=1, pos=0, dimRed=False):
    """
    Calculate the probability of measuring the value "measure".
    """
    p = projector_single_qubit_Zbasis(measure, N, pos)
    return np.linalg.norm((p * rho).tr())


def projector_single_qubit_Zbasis_dimRed(project, N=1, pos=0):
    # NOTE: This projector reduces the dimension of the state density matrix
    # is a rectangular matrix <x| not |x><x|
    if project != 0 and project != 1:
        raise ValueError("projector: measurement value invalid")

    p = qt.basis(2, project).dag()
    return tensor_single_operator(p, N, pos)


def projector_single_qubit_Zbasis(project, N=1, pos=0):
    # NOTE: This projector is |x><x|
    if project != 0 and project != 1:
        raise ValueError("projector: measurement value invalid")

    p = qt.basis(2, project)
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
