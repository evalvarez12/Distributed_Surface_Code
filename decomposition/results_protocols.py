"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import stabilizer
import protocols
import tools.names as names
import matplotlib.pyplot as plt
import noise_modeling
# Determine parameters
# NOTE: Realistic paramters
# ps = 0.006
# pm = 0.006
# pg = 0.006
# a0 = 83.33
# a1 = 1/3.
# eta = (0.1)*(0.03)*(0.8)
# theta = .63

# Improved parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 8
a1 = 1/30.
eta = 1/100
theta = .63

def env_error_rate(t, a):
    # Function to calculate the error to the enviroment for step of stabilizers
    # measurements
    x = a * t
    p_env = (1 + np.exp(-x))/2.
    return 1 - p_env

# Number of iterations for a average
iterations = 300
ignore_number = int(iterations/100)

# Initialize objects and define references
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
rho_ref2 = qt.bell_state('01') * qt.bell_state('01').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()
ghz3_ref = qt.ghz_state(3) * qt.ghz_state(3).dag()

# Stabilizer and error modeling stuff
stab_size = 4
parity = "X"

stab = stabilizer.Stabilizer(ps=ps, pm=pm, pg=pg)
model = noise_modeling.NoiseModel(stab_size, parity)
model.separate_basis_parity()

# Choi state for noise noise modeling
choi = model._choi_state_ket(stab_size)
choi = choi * choi.dag()
targets = list(range(stab_size))

# Lists to save results
FIDELITY = []
TIMES = []

# for eta in [1/30., 1/40., 1/50., 1/60., 1/70., 1/80.]:
# for a0 in [40., 30., 20., 10., 5., 2.]:
# for extra in [-20, -15, -10, -5, 0, 5, 10, 15, 20]:
# for a0 in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]:
# for s in [0]:
print("------> Var=", a0)
ghz = protocols.ghz4_single_simple(ps, pm, pg, eta, a0, a1, theta)
# Get average number of steps
fidelity = []
times = []
rho = ghz_ref*0
# check = collections.Counter({})
for i in range(iterations):
    print(i)
    # print(i)
    r, c = ghz.run()
    times += [c["time"]]

    r = stab.twirl_ghz(r)
    fidelity += [qt.fidelity(r, ghz_ref)]

    rho += r
rho = rho/iterations
times = np.array(times)
fidelity = np.array(fidelity)


# FIDELITY += [(np.average(fidelity), np.std(fidelity))]
# TIMES += [(np.average(times), np.std(times))]
# name = names.ghz(ps, pm, pg, eta, a0, a1, theta, 4, "DEMO")
# qt.qsave(rho, name)



indices_sorted = np.argsort(times)
t_sorted = times[indices_sorted]

if ignore_number != 0:
    t_max = t_sorted[:-ignore_number][-1]
    favg = np.average(fidelity[indices_sorted][:-ignore_number])
    fstd = np.std(fidelity[indices_sorted][:-ignore_number])
    tavg = np.average(times[indices_sorted][:-ignore_number])
    tstd = np.std(times[indices_sorted][:-ignore_number])
else:
    t_max = t_sorted[-1]
    favg = np.average(fidelity)
    fstd = np.std(fidelity)
    tavg = np.average(times)
    tstd = np.std(times)

print("F: ", favg, fstd)
print("T: ", tavg, tstd)
print("TIME_MAX:", t_max)


#############################################
#################### NOISE MODELING #########
p_res, rhos = stab.measure_ghz_stabilizer(choi, rho, targets, parity)
# Set channel output and make chi matrix
model.set_rho(rhos, p_res)
model.make_chi_matrix()

print("Total sum check: ", model.check_total_sum())
I_OK = model.chi["IIII_OK"]
I_NOK = model.chi["IIII_NOK"]
# The sum of all physical errors
E = (1 - model.chi["IIII_OK"] - model.chi["IIII_NOK"])/4.
print("Errors:")
print("OK: ", I_OK)
print("NOK: ", I_NOK)
print("E: ", E)
print("Env E: ", env_error_rate(t_max, a1)*4)
print("-------------")

############################################
############################################

# np.save("FIDELITY_DEMO_a0.npy", FIDELITY)
# np.save("TIME_DEMO_a0.npy", TIMES)


plt.plot(times[indices_sorted], fidelity[indices_sorted], 'bo')
if ignore_number != 0:
    vline = np.linspace(min(fidelity), max(fidelity))
    o = np.ones_like(vline)*(t_max)
    plt.plot(o, vline, 'r--')


plt.show()
