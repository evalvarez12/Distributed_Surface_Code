import matplotlib.pyplot as plt
import numpy as np

SIMPLE_fidelity = np.load("../decomposition/data/SIMPLE_fidelity.npy")
SIMPLE_steps = np.load("../decomposition/data/SIMPLE_steps.npy")
MEDIUM_fidelity = np.load("../decomposition/data/MEDIUM_fidelity.npy")
MEDIUM_steps = np.load("../decomposition/data/MEDIUM_steps.npy")
COMPLEX_fidelity = np.load("../decomposition/data/COMPLEX_fidelity.npy")
COMPLEX_steps = np.load("../decomposition/data/COMPLEX_steps.npy")
p_env = np.load("../decomposition/data/p_env.npy")

p_env = np.log(p_env)

plt.errorbar(p_env, SIMPLE_fidelity[:, 0], yerr=SIMPLE_fidelity[:, 1], fmt='o')
plt.errorbar(p_env, MEDIUM_fidelity[:, 0], yerr=MEDIUM_fidelity[:, 1], fmt='o')
plt.errorbar(p_env, COMPLEX_fidelity[:, 0], yerr=COMPLEX_fidelity[:, 1], fmt='o')

plt.figure()
# plt.plot(p_env, SIMPLE_steps[:, 0], 'o')
# plt.plot(p_env, MEDIUM_steps[:, 0],'o')
# plt.plot(p_env, COMPLEX_steps[:, 0], 'o')
plt.errorbar(p_env, SIMPLE_steps[:, 0], yerr=SIMPLE_steps[:, 1], fmt='o')
plt.errorbar(p_env, MEDIUM_steps[:, 0], yerr=MEDIUM_steps[:, 1], fmt='o')
plt.errorbar(p_env, COMPLEX_steps[:, 0], yerr=COMPLEX_steps[:, 1], fmt='o')
plt.show()
