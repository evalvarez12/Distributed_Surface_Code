"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import qutip as qt
import numpy as np
import protocols

# Determine parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 1
a1 = 1/80.
eta = 1/100
theta = .24

iterations = 10

# Initialize objects
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

F = []

print("-------------------PROTOCOL TEST------------------")
for i in range(iterations):
    circuit = protocols.purification_simple_4(ps, pm, pg, eta, a0, a1, theta)
    rho, operations = circuit.run(None)
    fidelity = qt.fidelity(rho, ghz_ref)
    print("F: ", fidelity)
    print("T: ", operations["time"])
    # print(operations)
    F += [fidelity]
    print("------------------------")

print("AVERAGE: ", np.average(F))
