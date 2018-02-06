"""
Functions to calculate the success probability of certain protocols.
Relies in 'projectors.py' to generate the required projectors.
"""
import projectors
import qutip as qt
import numpy as np


def single_sel(rho, N, ancillas_pos):
    P0 = (projector_single_qubit_Xbasis(0, N, ancillas_pos[0]) *
          projector_single_qubit_Xbasis(0, N, ancillas_pos[1]))
    P1 = (projector_single_qubit_Xbasis(1, N, ancillas_pos[0]) *
          projector_single_qubit_Xbasis(1, N, ancillas_pos[1]))

    P = P0 + P1
    p = (P * rho).tr().real
    return p

def epl(rho, N, ancillas_pos):
    P1 = (projector_single_qubit_Zbasis(1, N, ancillas_pos[0]) *
          projector_single_qubit_Zbasis(1, N, ancillas_pos[1]))

    p = (P1 * rho).tr().real
    return p

def double_sel(rho, N, ancillas_pos1, ancillas_pos2):
    P0a = (projector_single_qubit_Xbasis(0, N, ancillas_pos1[0]) *
           projector_single_qubit_Xbasis(0, N, ancillas_pos1[1]))
    P1a = (projector_single_qubit_Xbasis(1, N, ancillas_pos2[0]) *
           projector_single_qubit_Xbasis(1, N, ancillas_pos2[1]))

    P0b = (projector_single_qubit_Xbasis(0, N, ancillas_pos2[0]) *
           projector_single_qubit_Xbasis(0, N, ancillas_pos2[1]))
    P1b = (projector_single_qubit_Xbasis(1, N, ancillas_pos2[0]) *
           projector_single_qubit_Xbasis(1, N, ancillas_pos2[1]))

    # All the possible succes cases in the projectors
    P = P0a * P0b + P0a * P1b + P1a * P0b + P1a * P1b
    p =  (P * rho).tr().real
    return p


def ghz_4(rho, N, N_ancillas):
    """
    Requires that ancillas are on the last part of the entire state rho.
    """
    # Get a list of all even projectors for the 4 ancillas
    P_even = projectors.even_projectors(4)
    rest = N - N_ancillas
    P = sum(P_even)
    P = qt.tensor(qt.qeye([2]*rest), P)
    p =  (P * rho).tr().real
    return p


def ghz_correction_list(measurements, ghz_size):
    if ghz_size == 3:
        # Only two pairs of Bell states are used in this case
        if measurements == [0]:
            return qt.qeye([2]*3)
        elif measurements == [1]:
            return qt.rx(np.pi, 3, 2)

    if ghz_size == 4
        # There are two protocols in this case:
        #  1 or 2 ancilla bell pairs
        N_ancillas = len(measurements)
        if N_ancillas == 3:
            if sum(measurements) % 2 == 1:
                return qt.rx(np.pi, 4, 3) * qt.rx(np.pi, 4, 2)
            else:
                return qt.qeye([2]*4)
        elif N_ancillas == 4:
            if sum(measurements) % 2 == 1:
                return None
            no_correction_list = [0, 3, 12, 15]
            m_bin = int(''.join(map(str,measurements)),2)
            if m_bin in correction_list:
                return qt.qeye([2]*4)
            else:
                return qt.rx(np.pi, 4, 3) * qt.rx(np.pi, 4, 2)
