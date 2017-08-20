"""
Error models used to simulate the different types of decoherence effects.

author: Eduardo Villasenor
created-on: 05/06/17
"""

import qutip as qt
import numpy as np
import operations

sigmas = [qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()]

# TODO ask david if the order of the application noise gate is important

def single_qubit_gate_noise(rho, p, N=1, pos=0):
    res = rho*(1 - p)
    for i in sigmas[1:]:
        noise = operations.tensor_single_operator(i, N, pos)
        res += p/3.*(noise * rho * noise.dag())
    return res

def single_qubit_gate(rho, gate, p, N=1, pos=0):
    newRho = gate * rho * gate.dag()
    newRho = single_qubit_gate(newRho, p, N, pos)
    return newRho


def two_qubit_gate_noise(rho, p, N=2, pos1=0, pos2=1):
    res = rho*(1 - p)
    # Compute all the possible tensor products (i x j)
    for i in sigmas:
        for j in sigmas[1:]:
            noise = operations.tensor_operator([i, j], [pos1, pos2], N)
            res += p/15. * (noise * rho * noise.dag())

    # Do the missing last part (i x sig0)
    j = sigmas[0]
    for i in sigmas[1:]:
        noise = operations.tensor_operator([i, j], [pos1, pos2], N)
        res += p/15. * (noise * rho * noise.dag())
    return res

def two_qubit_gate(rho, gate, p, N=2, pos1=0, pos2=1):
    newRho = gate * rho * gate.dag()
    newRho = two_qubit_gate_noise(newRho, p, N, pos1, pos2)


def measurement_Zbasis(rho, p):
    res = measure
