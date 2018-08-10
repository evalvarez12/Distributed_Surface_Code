import matplotlib.pyplot as plt
import numpy as np

SIMPLE_fidelity = np.load("../decomposition/data/SIMPLE_fidelity.npy")
SIMPLE_steps = np.load("../decomposition/data/SIMPLE_steps.npy")
MEDIUM_fidelity = np.load("../decomposition/data/MEDIUM2_fidelity.npy")
MEDIUM_steps = np.load("../decomposition/data/MEDIUM2_steps.npy")
COMPLEX_fidelity = np.load("../decomposition/data/COMPLEX_fidelity.npy")
COMPLEX_steps = np.load("../decomposition/data/COMPLEX_steps.npy")
p_env = np.load("../decomposition/data/p_env.npy")


p_env = np.flip(np.logspace(-1, 2, len(p_env)), 0)

plt.errorbar(p_env, SIMPLE_fidelity[:, 0], yerr=SIMPLE_fidelity[:, 1], fmt='o', label="EPL")
plt.errorbar(p_env, MEDIUM_fidelity[:, 0], yerr=MEDIUM_fidelity[:, 1], fmt='o', label="SINGLE SELECTION")
plt.errorbar(p_env, COMPLEX_fidelity[:, 0], yerr=COMPLEX_fidelity[:, 1], fmt='o', label="DOUBLE SELECTION")

plt.ylabel(r"Fidelity", fontsize=13)
plt.xlabel(r"$a_0$", fontsize=13)
plt.xscale("log")
plt.legend(fontsize=13)
plt.figure()


plt.plot(p_env, SIMPLE_steps[:, 0], 'o', label="EPL")
plt.plot(p_env, MEDIUM_steps[:, 0],'o', label="SINGLE SELECTION")
plt.plot(p_env, COMPLEX_steps[:, 0], 'o', label="DOUBLE SELECTION")
# plt.errorbar(p_env, SIMPLE_steps[:, 0], yerr=SIMPLE_steps[:, 1], fmt='o')
# plt.errorbar(p_env, MEDIUM_steps[:, 0], yerr=MEDIUM_steps[:, 1], fmt='o')
# plt.errorbar(p_env, COMPLEX_steps[:, 0], yerr=COMPLEX_steps[:, 1], fmt='o')

plt.ylabel(r"Steps", fontsize=13)
plt.xlabel(r"$a_0$", fontsize=13)
plt.xscale("log")
plt.legend(fontsize=13)
plt.show()
