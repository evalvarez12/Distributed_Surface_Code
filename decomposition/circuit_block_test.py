"""
Test file for the functions on Circuit Block.

author: Eduardo Villasenor
created-on: 20/11/17
"""

import qutip as qt
import circuit_block
import numpy as np

# Determine parameters
# ps = 0.006
# pm = 0.006
# pg = 0.006
# a0 = 5.
# a1 = 1/80.
# eta = 1/100.
# theta = .24

# NOTE: Parameters as in Raja's thesis
ps = 0.006
pm = 0.006
pg = 0.006
a0 = 83.33
a1 = 1/3.
eta = (0.1)*(0.03)*(0.8)
theta = .24


# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
cb_ideal = circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4.)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()


print("------------------SINGLE SELECTION-------------------")
_, _, rho = cb.start_epl()
print("F initial: ", qt.fidelity(rho, rho_ref))
p, check, rho = cb.single_selection(rho, [0, 1], "Z")
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------DOUBLE SELECTION-------------------")
_, _, rho = cb.start_epl()
p, check, rho = cb.double_selection(rho, [0, 1], "X")
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------EPL protocol-------------------")
p, check, rho = cb.start_epl()
# p, check, rho = cb.double_selection(rho, [0, 1], "X")
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

# print("------------------NOISY BELL GENERATOR-------------------")
# p, check, rho = cb.start_bell_pair()
# print("p_success: ", p)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, rho_ref))
