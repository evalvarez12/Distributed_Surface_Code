"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import stabilizer
import matplotlib.pyplot as plt
import protocols
import circuit_block as circb
import circuit as circ

# Determine parameters
# NOTE: Realistic paramters
ps = 0.006
pm = 0.006
pg = 0.006
a0 = 83.33
a1 = 1/3.
eta = (0.1)*(0.03)*(0.8)
theta = 0.63

# Number of iterations for a average
iterations = 300
ignore_number = int(iterations/100)


cb = circb.Blocks(ps, pm, pg, eta, a0, a1, theta)

# First assemeble the small single selection circuit
single_sel = circ.Circuit(a0=a0, a1=a1,
                          circuit_block=cb.start_epl)
single_sel.add_circuit(circuit_block=cb.swap_pair,
                       pair=[0,1])
single_sel.add_circuit(circuit_block=cb.single_selection,
                       operation_qubits=[0, 1],
                       sigma="Z")





# Initialize objects and define references
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
rho_ref2 = qt.bell_state('01') * qt.bell_state('01').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()
ghz3_ref = qt.ghz_state(3) * qt.ghz_state(3).dag()
stab = stabilizer.Stabilizer(0, 0, 0)

# Lists to save results
FIDELITY = []
TIMES = []

for s in [0]:
    print("------> Var=", a0)

    circuit = protocols.pair_single_sel(ps, pm, pg, eta, a0, a1, theta)

    # Get average number of steps
    fidelity = []
    times = []
    fidelity2 = []
    times2 = []
    ps = []
    rho = rho_ref*0
    # check = collections.Counter({})
    for i in range(iterations):
        # print(i)
        r, c = circuit.runMC()

        r2, c2 = single_sel.run(None)

        times += [c["time"]]
        fidelity += [qt.fidelity(r, rho_ref)]

        times2 += [c2["time"]]
        fidelity2 += [qt.fidelity(r2, rho_ref)]

        rho += r
    rho = rho/iterations
    times = np.array(times)
    fidelity = np.array(fidelity)
    times2 = np.array(times2)
    fidelity2 = np.array(fidelity2)
    ps = np.array(ps)

    favg = np.average(fidelity)
    fstd = np.std(fidelity)
    tavg = np.average(times)
    tstd = np.std(times)

    f2avg = np.average(fidelity2)
    f2std = np.std(fidelity2)
    t2avg = np.average(times2)
    t2std = np.std(times2)


    print("F: ", favg, fstd)
    print("TIMES:", tavg, tstd)
    print("F2: ", f2avg, f2std)
    print("TIMES2:", t2avg, t2std)
    # FIDELITY += [(np.average(fidelity), np.std(fidelity))]
    # TIMES += [(np.average(times), np.std(times))]
    # name = names.ghz(ps, pm, pg, eta, a0, a1, theta, 4, "DEMO")
    # qt.qsave(rho, name)

# np.save("FIDELITY_DEMO_a0.npy", FIDELITY)
# np.save("TIME_DEMO_a0.npy", TIMES)
indices_sorted = np.argsort(times)

plt.plot(times[indices_sorted], 'bo')
plt.plot(fidelity[indices_sorted], 'ro')

indices_sorted2 = np.argsort(times2)

plt.plot(times2[indices_sorted2], 'bs')
plt.plot(fidelity2[indices_sorted2], 'rs')



plt.show()
