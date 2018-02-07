"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import protocols_tools
import protocols

# Determine parameters
# Very optimisitc gate and measurement paramterss
ps = 0.003
pm = 0.003
pg = 0.003

# Rajas optimistic entanglement and enviroment parameters
a1 = 1/3.
a0 = 83.33
eta = (0.1)*(0.03)*(0.8)

# My very optimistic parameters
# a0 = 8.0
# a1 = 1/80.
# eta = 1/100.

# Theta to optimize entanglement generation
theta = .24

# Number of iterations for a average
iterations = 50

# Initialize objects and define references
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
rho_ref2 = qt.bell_state('01') * qt.bell_state('01').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()
ghz3_ref = qt.ghz_state(3) * qt.ghz_state(3).dag()
pd = protocols_tools.Tools(0, 0, 0, 0)

# Lists to save results
FIDELITY = []
TIMES = []

# for eta in [1/30., 1/40., 1/50., 1/60., 1/70., 1/80.]:
# for a0 in [40., 30., 20., 10., 5., 2.]:
for extra in [-20, -15, -10, -5, 0, 5, 10, 15, 20]:
# for a0 in [10.0, 9.5, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0]:
    print("------> Var=", a0 + extra)
    ghz = protocols.BK_3(ps, pm, pg, eta, a0 + extra, a1, theta)
    # Get average number of steps
    fidelity = []
    times = []
    rho = ghz3_ref*0
    # check = collections.Counter({})
    for i in range(iterations):
        # print(i)
        _, c, r = ghz.run(None)
        times += [c["time"]]

        r = pd.twirl_ghz(r)
        fidelity += [qt.fidelity(r, ghz3_ref)]

        rho += r
    rho = rho/iterations
    print("F: ", np.average(fidelity), np.std(fidelity))
    print("TIMES:", np.average(times), np.std(times))
    FIDELITY += [(np.average(fidelity), np.std(fidelity))]
    TIMES += [(np.average(times), np.std(times))]
    # name = "ghz_3_eta_" + str(round(eta, 3))
    # name = names.ghz(ps, pm, pg, eta, a0, a1, theta, 3, "BK")
    # qt.qsave(rho, name)

np.save("Fidelity_Rajas_var=a0_bk_ghz3.npy", FIDELITY)
np.save("Times_Rajas_var=a0_bk_ghz3.npy", TIMES)
