"""
File to test the different entanglement generation protocols, when changing the
parameters.

author: Eduardo Villasenor
created-on: 04/01/18
"""

import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import protocols
import circuit_block

# Determine parameters
# # NOTE: Realistic paramters
ps = 0.006
pm = 0.006
pg = 0.006
a0 = 83.33
a1 = 1/3.
eta = (0.1)*(0.03)*(0.8)
theta = 0.63

# Improoved parameters
# ps = 0.006
# pm = 0.006
# pg = 0.006
# a0 = 5
# a1 = 1/80.
# eta = 1/100
# theta = .63


iterations = 30

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
a0s = np.logspace(0, 2, 8)

##################################################################
############### DATA IS OBTAINED HERE#############################

EPL = []
SINGLE = []
DOUBLE = []
for a0 in a0s:
    print(a0)
    cb.change_parameters(ps, pm, pg, eta, a0, a1, theta)

    f_epl = []
    f_single = []
    f_double = []
    for i in range(iterations):
        print(i)
        _, _, rho_epl = cb.start_epl()
        f_epl += [qt.fidelity(rho_epl, rho_ref)]

        _, _, rho_single = cb.single_selection(rho_epl, [0, 1], "Z")
        f_single += [qt.fidelity(rho_single, rho_ref)]

        _, _, rho_double= cb.double_selection(rho_epl, [0, 1], "Z")
        f_double += [qt.fidelity(rho_double, rho_ref)]

    EPL += [(np.average(f_epl), np.std(f_epl))]
    SINGLE += [(np.average(f_single), np.std(f_single))]
    DOUBLE += [(np.average(f_double), np.std(f_double))]

EPL = np.array(EPL)
SINGLE = np.array(SINGLE)
DOUBLE = np.array(DOUBLE)
# #
np.save("data/F_EPL.npy", EPL)
np.save("data/F_SINGLE.npy", SINGLE)
np.save("data/F_DOUBLE.npy", DOUBLE)


################### PLOT THE RESULTS ######################
# F = np.load("data/F_epl_params_raja.npy")
# T = np.load("data/T_epl_params_raja.npy")
# P = np.load("data/P_epl_params_raja.npy")

plt.figure()
plt.ylabel(r"$F$", fontsize=17)
plt.xlabel(r"$a_0$", fontsize=17)
plt.errorbar(a0s, EPL[:, 0], yerr=EPL[:, 1], fmt='o-', label=r"EPL")
plt.errorbar(a0s, SINGLE[:, 0], yerr=SINGLE[:, 1], fmt='o-', label=r"SINGLE")
plt.errorbar(a0s, DOUBLE[:, 0], yerr=DOUBLE[:, 1], fmt='o-', label=r"DOUBLE")
plt.legend(fontsize=17)
plt.xscale("log")
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)

plt.show()
