"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import stabilizer
import circuit_block
import circuit
import matplotlib.pyplot as plt

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
iterations = 1000
ignore_number = int(iterations/100)

# Initialize objects and define references
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
rho_ref2 = qt.bell_state('01') * qt.bell_state('01').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()
ghz3_ref = qt.ghz_state(3) * qt.ghz_state(3).dag()
stab = stabilizer.Stabilizer(0, 0, 0)


cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

# Lists to save results
FIDELITY = []
TIMES = []


# Initialize the circuit block object
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)

single_sel_simple = circuit.Circuit(a0=a0, a1=a1,
                                    circuit_block=cb.start_epl)

# single_sel_simple.add_circuit(circuit_block=cb.swap_pair,
#                               pair=[0, 1])
#
# single_sel_simple.add_circuit(circuit_block=cb.single_selection,
#                               operation_qubits=[0, 1],
#                               sigma="Z")


for s in [0]:
    print("------> Var=", a0)
    # Get average number of steps
    fidelity = []
    times = []
    ps = []
    rho = rho_ref*0
    # check = collections.Counter({})
    for i in range(iterations):
        # print(i)
        r, c = single_sel_simple.run(None)
        times += [c["time"]]
        fidelity += [qt.fidelity(r, rho_ref)]

        rho += r
    rho = rho/iterations
    times = np.array(times)
    fidelity = np.array(fidelity)
    ps = np.array(ps)

    favg = np.average(fidelity)
    fstd = np.std(fidelity)
    tavg = np.average(times)
    tstd = np.std(times)

    print("F: ", favg, fstd)
    print("TIMES:", tavg, tstd)
    # FIDELITY += [(np.average(fidelity), np.std(fidelity))]
    # TIMES += [(np.average(times), np.std(times))]
    # name = names.ghz(ps, pm, pg, eta, a0, a1, theta, 4, "DEMO")
    # qt.qsave(rho, name)

# np.save("FIDELITY_DEMO_a0.npy", FIDELITY)
# np.save("TIME_DEMO_a0.npy", TIMES)
indices_sorted = np.argsort(times)

plt.plot(times[indices_sorted], 'bo')
plt.plot(fidelity[indices_sorted], 'ro')



plt.show()
