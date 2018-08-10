"""
Functions to calculate the success probability of certain protocols.
Relies in 'projectors.py' and 'operations.py' to generate the required projectors.

author: Eduardo Villaseñor
created-on: 16/01/18
"""
import qutip as qt
import numpy as np
import faulttolerance.decomposition.tools.operations as ops
import faulttolerance.decomposition.tools.projectors as prj


def single_sel(rho, N, ancillas_pos):
    """
    Return the success probability of the single selection
    protocol for entanglement distillment.

    Parameters
    ----------
    rho : (densmat) density matrix.
    N : (int) system size.
    ancilla_pos : (int) position of qubit to be measured.
    """
    # Create the required projectors
    P0 = (ops.projector_single_qubit_Xbasis(0, N, ancillas_pos[0]) *
          ops.projector_single_qubit_Xbasis(0, N, ancillas_pos[1]))
    P1 = (ops.projector_single_qubit_Xbasis(1, N, ancillas_pos[0]) *
          ops.projector_single_qubit_Xbasis(1, N, ancillas_pos[1]))

    P = P0 + P1
    # Compute the probability
    p = (P * rho * P.dag()).tr()
    return p


def epl(rho, N, ancillas_pos):
    """
    Return the success probability of the single selection
    protocol for entanglement distillment.

    Parameters
    ----------
    rho : (densmat) density matrix.
    N : (int) system size.
    ancilla_pos : (int) position of qubit to be measured.
    """
    # Create the required projectors
    P = (ops.projector_single_qubit_Zbasis(1, N, ancillas_pos[0]) *
         ops.projector_single_qubit_Zbasis(1, N, ancillas_pos[1]))

    # Compute probability
    p = (P * rho * P).tr()
    return p


def double_sel(rho, N, ancillas_pos1, ancillas_pos2):
    """
    Return the success probability of the single selection
    protocol for entanglement distillment.

    Parameters
    ----------
    rho : (densmat) density matrix.
    N : (int) system size.
    ancilla_pos : (int) position of qubit to be measured.
    """
    # Create the required projectors
    P0a = (ops.projector_single_qubit_Xbasis(0, N, ancillas_pos1[0]) *
           ops.projector_single_qubit_Xbasis(0, N, ancillas_pos1[1]))
    P1a = (ops.projector_single_qubit_Xbasis(1, N, ancillas_pos1[0]) *
           ops.projector_single_qubit_Xbasis(1, N, ancillas_pos1[1]))

    P0b = (ops.projector_single_qubit_Xbasis(0, N, ancillas_pos2[0]) *
           ops.projector_single_qubit_Xbasis(0, N, ancillas_pos2[1]))
    P1b = (ops.projector_single_qubit_Xbasis(1, N, ancillas_pos2[0]) *
           ops.projector_single_qubit_Xbasis(1, N, ancillas_pos2[1]))

    # All the possible succes cases in the projectors
    P = (P0a * P0b) + (P0a * P1b) + (P1a * P0b) + (P1a * P1b)
    p = (P * rho * P.dag()).tr()
    return p


def ghz_4(rho):
    """
    Compute the probability of success for the special case when creating
    a GHZ state from 4 bell pairs.
    Requires that ancillas are on the last part of the entire state rho.

    Parameters
    ----------
    rho : (densmat) density matrix.
    """
    # Get a list of all even projectors for the 4 ancillas
    N = len(rho.dims[0])
    P_even = prj.even_projectors(4)
    rest = N - 4
    # Add an identity before the projector
    P = qt.tensor(qt.qeye([2]*rest), P_even)
    # Compute probability
    p = (P * rho * P.dag()).tr()
    return p

def ghz_3(rho):
    """
    Compute the probability of success for the special case when creating
    a GHZ state from 4 bell pairs.
    Requires that ancillas are on the last part of the entire state rho.

    Parameters
    ----------
    rho : (densmat) density matrix.
    """
    # Get a list of all even projectors for the 4 ancillas
    N = len(rho.dims[0])
    P_even = prj.even_projectors(3)
    rest = N - 3
    # Add an identity before the projector
    P = qt.tensor(qt.qeye([2]*rest), P_even)
    # Compute probability
    p = (P * rho * P.dag()).tr()
    return p
