import matplotlib.pyplot as plt
import numpy as np

SIMPLE_fidelity = np.load("../decomposition/data/SIMPLE_fidelity_BK.npy")
SIMPLE_steps = np.load("../decomposition/data/SIMPLE_steps_BK.npy")
MEDIUM_fidelity = np.load("../decomposition/data/MEDIUM_fidelity_BK.npy")
MEDIUM_steps = np.load("../decomposition/data/MEDIUM_steps_BK.npy")
COMPLEX_fidelity = np.load("../decomposition/data/COMPLEX_fidelity_BK.npy")
COMPLEX_steps = np.load("../decomposition/data/COMPLEX_steps_BK.npy")
p_env = np.load("../decomposition/data/p_env.npy")


plt.errorbar(p_env, SIMPLE_fidelity[:, 0], yerr=SIMPLE_fidelity[:, 1], fmt='o', label="SIMPLE")
plt.errorbar(p_env, MEDIUM_fidelity[:, 0], yerr=MEDIUM_fidelity[:, 1], fmt='o', label="MEDIUM")
plt.errorbar(p_env, COMPLEX_fidelity[:, 0], yerr=COMPLEX_fidelity[:, 1], fmt='o', label="COMPLEX")

plt.ylabel(r"Fidelity", fontsize=13)
plt.xlabel(r"$p_{env}$", fontsize=13)
plt.xscale("log")
plt.legend(fontsize=13)
plt.figure()


plt.plot(p_env, SIMPLE_steps[:, 0], 'o', label="SIMPLE")
plt.plot(p_env, MEDIUM_steps[:, 0],'o', label="MEDIUM")
plt.plot(p_env, COMPLEX_steps[:, 0], 'o', label="COMPLEX")
# plt.errorbar(p_env, SIMPLE_steps[:, 0], yerr=SIMPLE_steps[:, 1], fmt='o')
# plt.errorbar(p_env, MEDIUM_steps[:, 0], yerr=MEDIUM_steps[:, 1], fmt='o')
# plt.errorbar(p_env, COMPLEX_steps[:, 0], yerr=COMPLEX_steps[:, 1], fmt='o')

plt.ylabel(r"Steps", fontsize=13)
plt.xlabel(r"$p_{env}$", fontsize=13)
plt.xscale("log")
plt.legend(fontsize=13)
plt.show()
