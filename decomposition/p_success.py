"""
Functions to calculate the success probability of certain protocols.
Relies in 'projectors.py' and 'operations.py' to generate the required projectors.
"""
import projectors
import qutip as qt
import numpy as np
import operations as ops


def single_sel(rho, N, ancillas_pos):
    P_list = 2**ancillas_pos[0]


    P0 = (ops.projector_single_qubit_Xbasis(0, N, ancillas_pos[0]) *
          ops.projector_single_qubit_Xbasis(0, N, ancillas_pos[1]))
    P1 = (ops.projector_single_qubit_Xbasis(1, N, ancillas_pos[0]) *
          ops.projector_single_qubit_Xbasis(1, N, ancillas_pos[1]))

    P = P0 + P1
    p = (P * rho).tr().real
    return p

def epl(rho, N, ancillas_pos):
    P1 = (ops.projector_single_qubit_Zbasis(1, N, ancillas_pos[0]) *
          ops.projector_single_qubit_Zbasis(1, N, ancillas_pos[1]))

    p = (P1 * rho).tr().real
    return p

def double_sel(rho, N, ancillas_pos1, ancillas_pos2):
    P0a = (ops.projector_single_qubit_Xbasis(0, N, ancillas_pos1[0]) *
           ops.projector_single_qubit_Xbasis(0, N, ancillas_pos1[1]))
    P1a = (ops.projector_single_qubit_Xbasis(1, N, ancillas_pos2[0]) *
           ops.projector_single_qubit_Xbasis(1, N, ancillas_pos2[1]))

    P0b = (ops.projector_single_qubit_Xbasis(0, N, ancillas_pos2[0]) *
           ops.projector_single_qubit_Xbasis(0, N, ancillas_pos2[1]))
    P1b = (ops.projector_single_qubit_Xbasis(1, N, ancillas_pos2[0]) *
           ops.projector_single_qubit_Xbasis(1, N, ancillas_pos2[1]))

    # All the possible succes cases in the projectors
    P = P0a * P0b + P0a * P1b + P1a * P0b + P1a * P1b
    p =  (P * rho).tr().real
    return p


def ghz_4(rho):
    """
    Requires that ancillas are on the last part of the entire state rho.
    """
    # Get a list of all even projectors for the 4 ancillas
    N = len(rho.dims[0])
    P_even = projectors.even_projectors(4)
    rest = N - 4
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

    if ghz_size == 4:
        # There are two protocols in this case:
        #  1 or 2 ancilla bell pairs
        N_ancillas = len(measurements)
        if N_ancillas == 2:
            if sum(measurements) % 2 == 1:
                return qt.rx(np.pi, 4, 3) * qt.rx(np.pi, 4, 2)
            else:
                return qt.qeye([2]*4)
        elif N_ancillas == 4:
            if sum(measurements) % 2 == 1:
                return None
            no_correction_list = [0, 3, 12, 15]
            m_bin = int(''.join(map(str,measurements)),2)
            if m_bin in no_correction_list:
                return qt.qeye([2]*4)
            else:
                return qt.rx(np.pi, 4, 3) * qt.rx(np.pi, 4, 2)
