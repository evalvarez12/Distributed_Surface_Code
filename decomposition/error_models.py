"""
Error models used to simulate the decoherence effects due to imperfect operations.

author: Eduardo Villasenor
created-on: 05/06/17
"""

import qutip as qt
import numpy as np
import tools.operations as ops


def env_error_all(rho, a0, a1, t):
    """
    Apply environtmental dephasing on all qubits in the state rho.

    Parameters
    -----------
    rho : (densmat) the state to be dephased.
    a0 : (scalar) dephaising effect due to the entanglement generation attemts.
    a1 : (scalar) dephasing due to the environment interactions.
    t : (scalar) time in wich the dephasing ocurrs.
    """
    N = len(rho.dims[0])
    qubits = list(range(N))
    return env_error(rho, a0, a1, t, N, qubits)


def env_error(rho, a0, a1, t, N, qubits):
    """
    Apply environmental error on the selected qubits from the list.

    Parameters
    -----------
    rho : (densmat) the state to be dephased.
    a0 : (scalar) depolarizing effect due to the entanglement generation attemts.
    a1 : (scalar) depolarizing due to the environment interactions.
    t : (scalar) time in wich the dephasing ocurrs.
    N : (int) total size of the system
    qubits : (list) list with the qubits on which the error applies
    """
    for i in qubits:
        rho = env_error_single(rho, a0, a1, t, N, i)
    return rho


def env_error_single(rho, a0, a1, t, N=1, pos=0):
    """Apply environmental error on a single qubit from the state rho."""
    a = (a0 + a1)*t
    X = qt.rx(np.pi, N, pos)
    Y = qt.ry(np.pi, N, pos)
    Z = qt.rz(np.pi, N, pos)
    # ss = sigma_z * rho * sigma_z.dag()
    lamb = (1 + np.exp(-a * t))/4.
    rho = (lamb * rho + (1 - lamb) / 3. * (Z * rho * Z.dag() + X * rho * X.dag()
           + Y * rho * Y.dag()))
    return rho


def get_sigmas(N=1, pos=0):
    """Return set of single qubit X Y Z gates."""
    sigmas = [qt.qeye([2]*N), qt.rx(np.pi, N, pos), qt.ry(np.pi, N, pos),
              qt.rz(np.pi, N, pos)]
    return sigmas


def single_qubit_gate_noise(rho, p, N=1, pos=0):
    """Apply depolarizing noise due a single qubit operation."""
    res = rho*(1 - p)
    sigmas = get_sigmas(N, pos)
    for noise in sigmas[1:]:
        res += p/3.*(noise * rho * noise.dag())
    return res


def single_qubit_gate(rho, gate, p, N=1, pos=0):
    """
    Apply a single qubit with the corresponding noise.

    Parameters
    -----------
    rho : (densmat) the state to apply the noisy gate on
    gate : (operator) ideal gate to be applied
    p : (scalar) error rate
    N : (int) total system size
    pos : (int) position of the qubit in which the gate applies
    """
    new_rho = gate * rho * gate.dag()
    new_rho = single_qubit_gate_noise(new_rho, p, N, pos)
    return new_rho


def two_qubit_gate_noise(rho, p, N=2, pos1=0, pos2=1):
    """Depolarizing noise due to a noisy two qubit gate."""
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
    """
    Apply a two quit gate with its corresponding noise.

    Parameters
    -----------
    rho : (densmat) the state to apply the noisy gate on
    gate : (operator) ideal gate to be applied
    p : (scalar) error rate
    N : (int) total system size
    pos1 : (int) position of the first qubit in which the gate applies
    pos2 : (int) position of the second qubit in which the gate applies
    """
    rho = gate * rho * gate.dag()
    rho = two_qubit_gate_noise(rho, p, N, pos1, pos2)
    return rho


def measure_single_Zbasis_random(rho, p, N=1, pos=0):
    """
    Probabilisticly measure and collapse the state on the given position.

    Parameters
    -----------
    rho : (densmat) the state to apply the noisy gate on
    p : (scalar) measurement error rate
    N : (int) total system size
    pos : (int) position of the qubit to be measured
    """
    X = qt.rx(np.pi, N, pos)
    rho = (1 - p)*rho + p*(X * rho * X.dag())
    measurement, collapsed_rho = ops.random_measure_single_Zbasis(rho, N,
                                                                  pos, True)
    return measurement, collapsed_rho


def measure_single_Zbasis_forced(rho, p, project, N=1, pos=0):
    """
    Collapse the state on the given position and to the given projector.

    Parameters
    -----------
    rho : (densmat) the state to apply the noisy gate on
    p : (scalar) measurement error rate
    project : (scalar) 0 or 1, projector involved in the measurement collapse
    N : (int) total system size
    pos : (int) position of the qubit to be measured
    """
    X = qt.rx(np.pi, N, pos)
    rho = (1 - p)*rho + p*(X * rho * X.dag())
    p, collapsed_rho = ops.forced_measure_single_Zbasis(rho, N, pos,
                                                        project, True)
    return collapsed_rho


def measure_single_Xbasis_random(rho, p, N=1, pos=0):
    """Probabilisticly measure and collapse the state on the given position."""
    Z = qt.rz(np.pi, N, pos)
    rho = (1 - p)*rho + p*(Z * rho * Z.dag())
    measurement, collapsed_rho = ops.random_measure_single_Xbasis(rho, N,
                                                                  pos, True)
    return measurement, collapsed_rho


def measure_single_Xbasis_forced(rho, p, project, N=1, pos=0):
    """Collapse the state on the given position and to the given projector."""
    Z = qt.rz(np.pi, N, pos)
    rho = (1 - p)*rho + p*(Z * rho * Z.dag())
    p, collapsed_rho = ops.forced_measure_single_Xbasis(rho, N, pos,
                                                        project, True)
    return collapsed_rho


def bell_pair(p):
    """
    Return a noisy |phi+> Bell pair.

    Parameters
    -----------
    p : (scalar) error rate
    """
    a = qt.bell_state('00') * qt.bell_state('00').dag()
    b = qt.bell_state('01') * qt.bell_state('01').dag() \
        + qt.bell_state('10') * qt.bell_state('10').dag() \
        + qt.bell_state('11') * qt.bell_state('11').dag()
    W = (1-p)*a + p/3.*b
    # H = qt.snot(2, 1)
    # return H*W*H.dag()
    return W


def bell_pair_phi(p):
    """
    Return a noisy |phi+> Bell pair.

    Parameters
    -----------
    p : (scalar) error rate
    """
    a = qt.bell_state('10') * qt.bell_state('10').dag()
    b = qt.bell_state('01') * qt.bell_state('01').dag() \
        + qt.bell_state('00') * qt.bell_state('00').dag() \
        + qt.bell_state('11') * qt.bell_state('11').dag()
    W = (1-p)*a + p/3.*b
    # H = qt.snot(2, 1)
    # return H*W*H.dag()
    return W


def bell_pair_click(eta, theta):
    """
    Return a Bell pair generated using the single click protocol.
    Initial qubits are initialized in the states: sin(theta)|0> + cos(theta)|1>

    Parameters
    -----------
    eta : (scalar) detection efficiency
    theta : (scalar) determines how the initial qubits are initialized
    """

    s = np.sin(theta)**2
    r = ((1 - eta)*s)/(1 - eta*s)
    # TODO Check here
    state = qt.bell_state('10') * qt.bell_state('10').dag()
    noise = qt.tensor(qt.basis(2, 1), qt.basis(2, 1))
    noise = noise * noise.dag()
    return (1-r)*state + r*noise

def drift(state, pd):
    """
    Drift noise according to NN article.
    """
    N = len(state.dims[0])
    Z = qt.rz(np.pi, N, 0)
    state = (1-pd) * state + pd * Z * state * Z.dag()
    return state

def generate_noisy_ghz(F, N):
    """
    Generate a noisy GHZ state of a given fidelity.

    Parameters
    -----------
    F : (scalar) fidelity of the final state
    N : (int) size of the final GHZ
    """
    # p = 1 - F
    # nu = 4**N*p/(4**N - 1)
    ghz = qt.ghz_state(N)
    ghz = ghz * ghz.dag()
    rho = (F**2) * ghz + (1- F**2)/(2**N) * qt.qeye([2]*N)
    return rho
