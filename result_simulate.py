"""
File to test the surface code simulation.

created-on: 23/07/17
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt
import surface_code
import layers
import matching

def lambda_env(t, a0, a1):
    a = (a0 + a1)*t
    lamb = (1 + np.exp(-a * t))/2.
    return 1 - lamb


# Define the parameters
distance = 12
topology = "toric"
weights = [1, 1]

# Parameters for noisy measurement
ps = 0.008
pm = 0.008
pg = 0.008
eta = 0.0
a0 = 0.0
a1 = 0.0
protocol = "STRINGENT"
theta = .0
PERFECT_LAST_ROUND = False

p = 0.12
q = 0.12
iterations = 1
cycles = 1

# Initialize objects
fail_rate = 0
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)
# sc.init_error_obj(topology, ps, pm, pg, eta, a0, a1, theta, protocol)

# Choose a measurement protocol
# sc.select_measurement_protocol(0, 0, "local")

# Perform measurements
for i in range(iterations):

    # Errors and measurements
    # Random errors
    if q != 0:
        for t in range(cycles):
            sc.apply_qubit_error(p, 0)
            sc.measure_all_stabilizers()
            sc._stabilizer_lie("S", q)
            lc.add()
        sc.plot("star")
        plt.savefig('sc.pdf', format='pdf', dpi=300)
        sc.measure_all_stabilizers()
        lc.add()
    else:
        sc.apply_qubit_error(p, 0)
        sc.measure_all_stablizers()
        lc.add()

    # Noisy measurements
    # for t in range(cycles):
    #     sc.noisy_measurement_cycle()
    #     lc.add()
    # sc.plot_all()
    # qubits_copy = sc.qubits.copy()
    #
    # sc.measure_all_stabilizers()
    # lc.add()

    # Decode
    lc.decode()


    # Round of perfect detection to eliminate stray errors
    if PERFECT_LAST_ROUND:
        lc.reset()
        sc.measure_all_stablizers()
        lc.add()
        anyons_star, anyons_plaq = lc.find_anyons_all()
        match_star = matching.match(distance, anyons_star, topology,
                                    "star", time=0, weights=weights)
        match_plaq = matching.match(distance, anyons_plaq, topology,
                                    "plaq", time=0, weights=weights)
        sc.correct_error("star", match_star, cycles)
        sc.correct_error("plaq", match_plaq, cycles)

    # Check for errors in decoding and correcting
    sc.measure_all_stabilizers()
    if (sc.qubits[:, sc.tags != "Q"] == -1).any():
        print("FAILURE CORRECTING")

    # Measure logical qubit
    logical = sc.measure_logical()

    # sc.plot_all()
    sc.plot("star")
    # plt.show()
    plt.savefig('sc_corrected.pdf', format='pdf', dpi=300)

    # Code to check when a logical error happens
    print(logical)
    if -1 in logical[0] or -1 in logical[1]:
        fail_rate += 1
        # sc.plot("star")
        # sc.plot("plaq")


        # pre_star = pre_correction[0].copy()
        # pre_star[sc.tags == "S"] *= 2
        # pre_star[sc.tags == "P"] = 1
        # pre_plaq = pre_correction[0].copy()
        # pre_plaq[sc.tags == "Q"] = pre_correction[1, sc.tags == "Q"]
        # pre_plaq[sc.tags == "P"] *= 2
        # pre_plaq[sc.tags == "S"] = 1

        # Return data to plot
        # return data, self.cmap, self.cmap_norm
        # plt.figure()
        # plt.imshow(pre_star)
        # plt.figure()
        # plt.imshow(pre_plaq)
        # plt.colorbar()
        # plt.show()

        # plt.show()

    lc.reset()
    sc.reset()

fail_rate = fail_rate/float(iterations)
print("FAIL RATE: ", fail_rate)



# print(sc.errors.n_NOK)
# print(sc.errors.n_E)

plt.show()


#
# x = sc.qubits.copy()
# x.fill(1)
# s = sc.stars.copy()
# s1 = sc.stars_round1.copy()
# s2 = sc.stars_round2.copy()
# p1 = sc.plaqs_round1.copy()
# p2 = sc.plaqs_round2.copy()
#
# sq1 = sc._stabilizer_qubits_bulk(s1)
# sq = sc._stabilizer_qubits_bulk(s)
#
# e = np.ones_like(sq).transpose((1, 0, 2))*3
# x[:, sq[:, 0], sq[:, 1]] *= e
