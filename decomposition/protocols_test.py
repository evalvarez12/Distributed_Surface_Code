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
a0 = 5
a1 = 1/80.
eta = 1/100
theta = .63

# NOTE: Parameters as in Raja's thesis
# ps = 0.006
# pm = 0.006
# pg = 0.006
# a0 = 83.33
# a1 = 1/3.
# eta = (0.1)*(0.03)*(0.8)
# theta = .24

iterations = 30

# Initialize objects
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

F = []
T = []
print("-------------------PROTOCOL TEST------------------")
for i in range(iterations):
    
    circuit = protocols.pair_single_sel(ps, pm, pg, eta, a0, a1, theta)
    rho, operations = circuit.run(None)
    fidelity = qt.fidelity(rho, ghz_ref)
    print("F: ", fidelity)
    print("T: ", operations["time"])
    # print(operations)
    # print(operations)
    F += [fidelity]
    T += [operations["time"]]
    # print("------------------------")

print("AVERAGE: ", np.average(F))
print("TIME: ", np.average(T))
