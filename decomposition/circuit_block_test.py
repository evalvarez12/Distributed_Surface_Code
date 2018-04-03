"""
Test file for the functions on Circuit Block.

author: Eduardo Villasenor
created-on: 20/11/17
"""

import qutip as qt
import circuit_block
import numpy as np
import error_models as errs


def p_ref(f): return (f**2 +2*f*(1-f)/3 + 5*((1-f)/3)**2)
def single_selection(F1, F2):
    F = F1*F2 + (1 - F1)*(1 - F2)/9
    F = F/(F1*F2 + F1*(1 - F2)/3 + F2*(1 - F1)/3 + 5*(1 - F1)*(1 - F2)/9)
    return F

# Determine parameters
ps = 0.006
pm = 0.006
pg = 0.006
a0 = 5.
a1 = 1/80.
eta = 1/100.
theta = .24

# NOTE: Parameters as in Raja's thesis
# ps = 0.006
# pm = 0.006
# pg = 0.006
# a0 = 83.33
# a1 = 1/3.
# eta = (0.1)*(0.03)*(0.8)
# theta = .24


# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
cb_ideal = circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4.)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()

print("------------------TEST SWAP--------------------------")
_, _, rho = cb.start_epl()
print("F initial: ", qt.fidelity(rho, rho_ref))
rho = cb._swap_pair(rho, [0, 1])
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------SINGLE SELECTION-------------------")
rho = errs.bell_pair(0.1)
print("F initial: ", qt.fidelity(rho, rho_ref)**2)
rho = qt.tensor(rho, rho)
p, check, rho = cb_ideal.single_selection_ops(rho, [0, 1], [2, 3],  "Z")
print("p_success: ", p, p_ref(0.9))
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref)**2, single_selection(.9, .9))

print("------------------DOUBLE SELECTION-------------------")
rho = errs.bell_pair(0.1)
print("F initial: ", qt.fidelity(rho, rho_ref)**2)
rho = qt.tensor(rho, rho, rho)
p, check, rho = cb_ideal.double_selection_ops(rho, [0, 1], [2, 3], [4, 5], "Z")
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref)**2)

# print("------------------DOUBLE SELECTION2-------------------")
# p, check, rho = cb.double_selection22("Z")
# print("p_success: ", p)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, rho_ref))
#
#
# print("------------------EPL protocol-------------------")
# p, check, rho = cb.start_epl()
# # p, check, rho = cb.double_selection(rho, [0, 1], "X")
# print("p_success: ", p)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, rho_ref))

# print("------------------NOISY BELL GENERATOR-------------------")
# p, check, rho = cb.start_bell_pair()
# print("p_success: ", p)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, rho_ref))
