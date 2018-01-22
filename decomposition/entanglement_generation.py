"""
File to test the different entanglement generation protocols, when changing the
parameters.

author: Eduardo Villasenor
created-on: 04/01/18
"""

import qutip as qt
import circuit_block
import circuit
import numpy as np
import matplotlib.pyplot as plt


# Determine parameters
# NOTE: Parameters as in Raja's thesis
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 83.33
a1 = 1/3.
eta = (0.1)*(0.03)*(0.8)
theta = np.pi/4.

iterations = 200

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
cb_ideal = circuit_block.Blocks(0, 0, 0, 1, 0, 0, np.pi/4.)
epl = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_epl)

rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
thetas = np.linspace(0.1, np.pi/2. - 0.1, 20)

##################################################################
############### DATA IS OBTAINED HERE#############################
#
# F = []
# T = []
# P = []
# for theta in thetas:
#     print(theta)
#     cb.change_parameters(ps, pm, pg, eta, a0, a1, theta)
#     f = []
#     t = []
#     p = []
#     for i in range(iterations):
#         # p_success, check, rho = epl.run(None)
#         p_success, check, rho = cb.start_epl(None)
#
#         p += [p_success]
#         f += [qt.fidelity(rho, rho_ref)]
#         t += [check["time"]]
#     T += [(np.average(t), np.std(t))]
#     F += [(np.average(f), np.std(f))]
#     P += [(np.average(p), np.std(p))]
#
# F = np.array(F)
# T = np.array(T)
# P = np.array(P)
# #
# np.save("F_good.npy", F)
# np.save("T_good.npy", T)
# np.save("P_good.npy", P)

#################################################################
#################################################################


# Define functions that return the theoretical values for reference

# One click protocol success probability
def p_one_click(t):
    c = np.cos(t)**2
    s = np.sin(t)**2
    return 2*s*(c*eta + s*eta*(1-eta))

# One click protocol fidelity
def f_one_click(t):
    c = np.cos(t)**2
    s = np.sin(t)**2
    return c/(c + s*(1-eta))

# Barret-Kok protocol success probability
def p_bk(t):
    return (f_one_click(t)**2)*eta**2

# Single selection resulting fidelity of purification with two states with
# fidelities F1 and F2
def single_selection(F1, F2):
    F = F1*F2 + (1 - F1)*(1 - F2)/9
    F = F/(F1*F2 + F1*(1 - F2)/3 + F2*(1 - F1)/3 + 5*(1 - F1)*(1 - F2)/9)
    return F

################### PLOT THE RESULTS ######################
F = np.load("data/F_epl_params_raja.npy")
T = np.load("data/T_epl_params_raja.npy")
P = np.load("data/P_epl_params_raja.npy")
#
plt.figure()
plt.ylabel(r"$F$", fontsize=17)
plt.xlabel(r"$\theta$", fontsize=17)
plt.errorbar(thetas, F[:, 0], yerr=F[:, 1], fmt='go-', label=r"EPL")
plt.plot(thetas, f_one_click(thetas), 'r-',label="One click")
plt.legend(fontsize=17)
plt.xticks(fontsize=17)


plt.figure()
plt.ylabel(r"$t(s)$", fontsize=17)
plt.plot(thetas, T[:, 0]/P[:,0], 'go-', label=r"EPL")
plt.plot(thetas, T[:, 0]/2., 'r-',label="One click")
plt.xlabel(r"$\theta$", fontsize=17)
# plt.legend(fontsize=17)
plt.xticks(fontsize=17)

#
plt.figure()
plt.ylabel(r"$p_{success}$", fontsize=17)
plt.xlabel(r"$\theta$", fontsize=17)
plt.errorbar(thetas, P[:, 0], yerr=P[:, 1], fmt='go-', label=r"EPL")
plt.plot(thetas, 100*p_one_click(thetas), 'r-', label="EOne click")
plt.plot(thetas, 1000*p_one_click(thetas)*P[:, 0], 'k--' )
plt.plot(thetas, 10000*p_bk(thetas), 'k--' )
plt.xticks(fontsize=17)

# plt.legend()
# plt.plot(thetas, 500000*p_bk(thetas))


plt.show()
