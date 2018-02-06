"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import protocols

# Determine parameters
ps = 0.006
pm = 0.006
pg = 0.006
a0 = 0
a1 = 0.
eta = 1
theta = .24

iterations = 1

# Initialize objects
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()



'''PROTOCOLS START HERE'''
print("-------------------PROTOCOL TEST------------------")
circuit = protocols.EPL_4(ps, pm, pg, eta, a0, a1, theta)
p, ops, rho = circuit.run(None)
fidelity = qt.fidelity(rho, ghz_ref)
print(fidelity)
