"""
Functions to apply the required corrections when generating GHZ states under
the different numbers of Bell pairs and GHZ sizes.

author: Eduardo Villasenor
created-on: 02/05/18
"""
import qutip as qt
import numpy as np


def correction(measurements, ghz_size, N, operation_pos):
    """
    Correction operations necesary given the result of the measurements and
    the GHZ state size.

    Parameters
    ----------
    measuements : (int) Measurements results of the collaped ancillas ex. [0, 1, 1, 0].
    N : (int) system size.
    operation_pos : (list) position of qubit(s) in which the corrections apply.
    """
    if ghz_size == 4:
        return correction_ghz4(measurements, N, operation_pos)
    elif ghz_size == 3:
        return correction_ghz3(measurements, N, operation_pos[0])

def correction_ghz3(measurements, N, operation_pos):
    """
    Correction opertions for the GHZ state of size 3.
    On the case when 3 bell pairs are collapsed into the GHZ.
    """

    # Only one qubit in the Bell states are used is measured
    # Result 0 requires no correction, result 1 requires to apply a Pauli X
    if measurements == [0]:
        return []
    elif measurements == [1]:
        return [qt.rx(np.pi, N, operation_pos)]


def correction_ghz4(measurements, N, operation_pos):
    # print("Measurements: ", measurements)
    # There are two protocols in this case:
    #  1 or 2 ancilla bell pairs
    N_ancillas = len(measurements)
    # If 1 bell pair is used: Then even measument sum requires no correction.
    # Odd measurement requires X correction
    if N_ancillas == 2:
        if sum(measurements) % 2 == 1:
            return [qt.rx(np.pi, N, operation_pos[0]),
                    qt.rx(np.pi, N, operation_pos[1])]
        else:
            return []
    # If two pairs are used then odd measurement means the state collapsed into
    # something mixed, as a result of entanglement disitillation.
    # Even result may o no require correction, do calculation by hand
    elif N_ancillas == 4:
        if sum(measurements) % 2 == 1:
            return None

        # Results that require no correction.
        no_correction_list = [0, 3, 12, 15]
        # Transform measurement result into int by binary transform
        m_bin = int(''.join(map(str, measurements)), 2)
        if m_bin in no_correction_list:
            return []
        else:
            return [qt.rx(np.pi, N, operation_pos[0]),
                    qt.rx(np.pi, N, operation_pos[1])]
