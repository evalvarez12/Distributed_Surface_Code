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
distance = 10
topology = "toric"
weights = [1, 1]

# Parameters for noisy measurement
ps = 0.003
pm = 0.003
pg = 0.003
eta = 0.01
a0 = 12.0
a1 = 1/80.
protocol = "GHZ"
theta = .24
NOISY_MEASUREMENT = True

p = 0.01
q = 0.01
iterations = 1
cycles = 5

# Initialize objects
fail_rate = 0
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)
sc.init_error_obj(topology, ps, pm, pg, eta, a0, a1, theta, protocol)

# Set time for each GHZ generation
time = 0.30347
lamb = lambda_env(time, 0, a1)
print("LAMBDA: ", lamb)
# Perform measurements
for i in range(iterations):

    # Errors and measurements
    # Random errors
    if q != 0:
        for t in range(cycles):
            sc.apply_qubit_error(p, 0)
            sc.measure_all_stablizers()
            sc.apply_measurement_error(q)
            lc.add()
    else:
        sc.apply_qubit_error(p, 0)
        sc.measure_all_stablizers()
        lc.add()

    # Noisy measurements
    # for t in range(cycles):
    #     # sc.noisy_measurement("star")
    #     # sc.noisy_measurement("plaq")
    #     sc.noisy_measurement_cycle(lamb)
    #     lc.add()

    # Get anyons
    anyons_star, anyons_plaq = lc.find_anyons_all()

    # Decode
    match_star = matching.match(distance, anyons_star, topology,
                                "star", time=cycles, weights=[1, 1])
    match_plaq = matching.match(distance, anyons_plaq, topology,
                                "plaq", time=cycles, weights=[1, 1])
    sc.plot("star")
    sc.plot("plaq")
    pre_correction = sc.qubits.copy()

    # Apply corrections
    sc.correct_error("star", match_star)
    sc.correct_error("plaq", match_plaq)

    # Round of perfect detection to eliminate stray errors
    if q!= 0 or NOISY_MEASUREMENT:
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
    sc.measure_stabilizer_type("star")
    sc.measure_stabilizer_type("plaq")
    if (sc.qubits[:, sc.tags != "Q"] == -1).any():
        print("FAILURE CORRECTING")
        fail_rate = -9999

    # Measure logical qubit
    logical = sc.measure_logical()

    sc.plot("star")
    sc.plot("plaq")
    plt.show()

    # Code to check when a logical error happens
    if -1 in logical[0] or -1 in logical[1]:
        fail_rate += 1
        print(logical)
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
print("FAILE RATE: ", fail_rate)

plt.show()
