"""
Functions to calculate the success probability of certain protocols.
Relies in 'projectors.py' and 'operations.py' to generate the required projectors.
"""
import qutip as qt
import numpy as np


def correction(measurements, ghz_size, N, ancillas_pos):
    if ghz_size == 4:
        return correction_ghz4(measurements, N, ancillas_pos)
    elif ghz_size == 3:
        return correction_ghz3(measurements, N, ancillas_pos)

def correction_ghz3(measurements, N, ancilla_pos):
    # Only two pairs of Bell states are used in this case
    if measurements == [0]:
        return []
    elif measurements == [1]:
        return [qt.rx(np.pi, N, ancilla_pos)]


def correction_ghz4(measurements, N, ancillas_pos):
        # There are two protocols in this case:
        #  1 or 2 ancilla bell pairs
        N_ancillas = len(measurements)
        if N_ancillas == 2:
            if sum(measurements) % 2 == 1:
                return [qt.rx(np.pi, N, ancillas_pos[0]),
                        qt.rx(np.pi, N, ancillas_pos[1])]
            else:
                return []
        elif N_ancillas == 4:
            if sum(measurements) % 2 == 1:
                return None

            no_correction_list = [0, 3, 12, 15]
            m_bin = int(''.join(map(str, measurements)), 2)
            if m_bin in no_correction_list:
                return []
            else:
                return [qt.rx(np.pi, N, ancillas_pos[0]),
                        qt.rx(np.pi, N, ancillas_pos[1])]
