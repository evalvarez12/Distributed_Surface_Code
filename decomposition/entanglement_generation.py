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

# print("------------------EPL protocol-------------------")
# p, check, rho = cb.start_epl()
# # p, check, rho = cb.double_selection(rho, [0, 1], "X")
# print("p_success: ", p)
# print("check: ", check)
# print("F: ", qt.fidelity(rho, rho_ref))

F = []
T = []
P = []
thetas = np.linspace(0.1, np.pi/2. - 0.1, 30)
for theta in thetas:
    print(theta)
    cb.change_parameters(ps, pm, pg, pn, pd, a0, a1, theta)
    f = []
    t = []
    p = []
    for i in range(iterations):
        # p_success, check, rho = epl.run(None)
        p_success, check, rho = cb.start_epl(None)

        p += [p_success]
        f += [qt.fidelity(rho, rho_ref)]
        t += [check["time"]]
    T += [(np.average(t), np.std(t))]
    F += [(np.average(f), np.std(f))]
    P += [(np.average(p), np.std(p))]

F = np.array(F)
T = np.array(T)
P = np.array(P)

np.save("F.npy", F)
np.save("T.npy", T)
np.save("P.npy", P)


plt.figure()
plt.ylabel(r"$F$")
plt.errorbar(thetas, F[:, 0], yerr=F[:, 1], fmt='bo-')

plt.figure()
plt.ylabel(r"$t$")
plt.plot(thetas, T[:, 0], 'ro-')

plt.figure()
plt.ylabel(r"$p_success$")
plt.errorbar(thetas, P[:, 0], yerr=P[:, 1], fmt='go-')

plt.show()
