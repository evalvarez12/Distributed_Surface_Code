import qutip as qt
import circuit_block
import circuit
import error_models as errs
import operations as ops
import numpy as np
import matplotlib.pyplot as plt
# Determine parameters
ps = 0.006
pm = 0.006
pg = 0.006
pn = 0.0
a0 = 83
a1 = 1/3.
pd = 1/2000.
theta = np.pi/4.

iterations = 200

# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, pn, pd, a0, a1, theta)
cb_ideal = circuit_block.Blocks(0, 0, 0, 0, 1, 0, 0, np.pi/4.)
epl = circuit.Circuit(a0=a0, a1=a1, circuit_block=cb.start_epl)

rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()

##################################################################
##################################################################

# F = []
# T = []
# P = []
thetas = np.linspace(0.1, np.pi/2. - 0.1, 30)
# for theta in thetas:
#     print(theta)
#     cb.change_parameters(ps, pm, pg, pn, pd, a0, a1, theta)
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
#
# np.save("F.npy", F)
# np.save("T.npy", T)
# np.save("P.npy", P)

#################################################################
#################################################################

# One click
def p_one_click(t):
    c = np.cos(t)**2
    s = np.sin(t)**2
    return 2*s*(c*pd + s*pd*(1-pd))

def f_one_click(t):
    c = np.cos(t)**2
    s = np.sin(t)**2
    return c/(c + s*(1-pd))

def p_bk(t):
    return (f_one_click(t)**2)*pd**2

F = np.load("F2.npy")
T = np.load("T2.npy")
P = np.load("P2.npy")

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
plt.plot(thetas, 500*p_one_click(thetas), 'r-', label="EOne click")
plt.plot(thetas, 1000*p_one_click(thetas)*P[:, 0], 'k--' )
plt.xticks(fontsize=17)

# plt.legend()
# plt.plot(thetas, 500000*p_bk(thetas))


plt.show()
