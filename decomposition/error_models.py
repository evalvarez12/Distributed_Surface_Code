"""
Error models used to simulate the different types of decoherence effects.

author: Eduardo Villasenor
created-on: 05/06/17
"""

import qutip as qt
import numpy as np
# from . import operations as ops
import operations as ops

# sigmas = [qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()]

# TODO ask david if the order of the application noise gate is important

def env_dephasing_all(rho, n_steps, e_spin):
    N = len(rho.dims[0])
    qs = list(range(N))
    return env_dephasing(rho, n_steps,e_spin, N, qs)

def env_dephasing(rho, n_steps, e_spin, N, qubits):
    for i in qubits:
        rho = env_dephasing_single(rho, n_steps, e_spin, N, i)
    return rho


def env_dephasing_single(rho, n_steps, e_spin, N=1, pos=0):
    # e_spin = TRUE or FALSE if electron spin is being used or no
    # NOTE: Time it takes to make every operation 2 micro seg
    # t_step = 2e-6
    # a0 = 1/2000
    # a1 = 1/3
    a = (5e-4) * e_spin + 1/3 * (2e-6)
    sigma_z = qt.rz(np.pi, N, pos)
    ss = sigma_z * rho * sigma_z.dag()
    lamb = np.exp(-a * n_steps)
    return (1 + lamb)/2 * rho + (1 - lamb)/2 * sigma_z * rho * sigma_z.dag()


def get_sigmas(N=1, pos=0):
    sigmas = [qt.qeye([2]*N), qt.rx(np.pi, N, pos), qt.ry(np.pi, N, pos),
              qt.rz(np.pi, N, pos)]
    return sigmas


def single_qubit_gate_noise(rho, p, N=1, pos=0):
    res = rho*(1 - p)
    sigmas = get_sigmas(N, pos)
    for noise in sigmas[1:]:
        res += p/3.*(noise * rho * noise.dag())
    return res


def single_qubit_gate(rho, gate, p, N=1, pos=0):
    new_rho = gate * rho * gate.dag()
    new_rho = single_qubit_gate_noise(new_rho, p, N, pos)
    return new_rho


def two_qubit_gate_noise(rho, p, N=2, pos1=0, pos2=1):
    # Remove the extra added p/15 due to
    # the 2 identities in the sigmas list
    res = rho*(1 - p)
    sigmas_pos1 = get_sigmas(N, pos1)
    sigmas_pos2 = get_sigmas(N, pos2)

    # Compute all the possible tensor products (i x j)
    for i in sigmas_pos1:
        for j in sigmas_pos2[1:]:
            noise = i * j
            res += p/15. * (noise * rho * noise.dag())

    # Do the missing last part (i x sig0)
    j = sigmas_pos2[0]
    for i in sigmas_pos1[1:]:
        noise = i * j
        res += p/15. * (noise * rho * noise.dag())
    return res


def two_qubit_gate(rho, gate, p, N=2, pos1=0, pos2=1):
    rho = gate * rho * gate.dag()
    rho = two_qubit_gate_noise(rho, p, N, pos1, pos2)
    return rho


def measure_single_Zbasis_random(rho, p, N=1, pos=0):
    X = qt.rx(np.pi, N, pos)
    rho = (1 - p)*rho + p*(X * rho * X.dag())
    measurement, collapsed_rho = ops.random_measure_single_Zbasis(rho, N,
                                                                  pos, True)
    return measurement, collapsed_rho


def measure_single_Zbasis_forced(rho, p, project, N=1, pos=0):
    X = qt.rx(np.pi, N, pos)
    rho = (1 - p)*rho + p*(X * rho * X.dag())
    p, collapsed_rho = ops.forced_measure_single_Zbasis(rho, N, pos, True)
    return p, collapsed_rho


def measure_single_Xbasis_random(rho, p, N=1, pos=0):
    Z = qt.rz(np.pi, N, pos)
    rho = (1 - p)*rho + p*(Z * rho * Z.dag())
    measurement, collapsed_rho = ops.random_measure_single_Xbasis(rho, N,
                                                                  pos, True)
    return measurement, collapsed_rho


def measure_single_Xbasis_forced(rho, p, project, N=1, pos=0):
    Z = qt.rz(np.pi, N, pos)
    rho = (1 - p)*rho + p*(Z * rho * Z.dag())
    p, collapsed_rho = ops.forced_measure_single_Xbasis(rho, N, pos,
                                                        project, True)
    return p, collapsed_rho


def bell_pair(p):
    a = qt.bell_state('00') * qt.bell_state('00').dag()
    b = qt.bell_state('01') * qt.bell_state('01').dag() \
        + qt.bell_state('10') * qt.bell_state('10').dag() \
        + qt.bell_state('11') * qt.bell_state('11').dag()
    W = (1-p)*a + p/3.*b
    # H = qt.snot(2, 1)
    # return H*W*H.dag()
    return W
