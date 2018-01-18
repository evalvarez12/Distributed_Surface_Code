import numpy as np
import qutip as qt
import entanglement_generation as egen


# Determine parameters
p_em = 0.1
p_ps = 0.03
p_det = 0.8
theta = np.pi/4.

# Initialize object
eg = egen.EntanglementGeneration(p_em=p_em, p_ps=p_ps, p_det=p_det,
                                 theta=theta)

# Reference state
rho_ref = qt.bell_state('01')

time, state = eg.generate_bell_pair_BK()
print("Time: ", time)
print("fidelity: ", qt.fidelity(state, rho_ref))
