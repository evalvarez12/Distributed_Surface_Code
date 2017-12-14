import qutip as qt
import circuit_block
import circuit
import error_models as errs
import operations as ops
import numpy as np
import matplotlib.pyplot as plt
# Determine parameters
ps = 0.005
pm = 0.005
pg = 0.005
pn = 0.1
a0 = 20
a1 = 1/3.
pd = 1/1000.
theta = np.pi/4.

iterations = 10

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
thetas = np.linspace(0.1, np.pi/2., 10)
for theta in thetas:
    print(theta)
    cb.change_parameters(ps, pm, pg, pn, pd, a0, a1, theta)
    f = []
    t = []
    for i in range(iterations):
        p, check, rho = epl.run(None)

        f += [qt.fidelity(rho, rho_ref)]
        t += [check["time"]]
    T += [(np.average(t), np.std(t))]
    F += [(np.average(f), np.std(f))]

F = np.array(F)
T = np.array(T)

plt.figure()
plt.errorbar(thetas, F[:, 0], yerr=F[:, 1], fmt='o-')

plt.figure()
plt.plot(thetas, T[:, 0], 'o-')

plt.show()
