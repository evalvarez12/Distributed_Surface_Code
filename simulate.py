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

# Define the parameters
distance = 10
topology = "planar"
time_steps = 20
weights = [1, 1]

# Parameters for noisy measurement
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.0
protocol = "LOCAL"

# Initialize objects
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)
sc.init_error_obj(ps, pm, pg, pn, protocol)


# Perform measurements
for i in range(time_steps):
    sc.apply_qubit_error(.03, .0)
    sc.measure_all_stablizers()
    # sc._stabilizer_lie("S", .00)
    # sc.noisy_measurement("star")
    lc.add()


    time = lc.get_time()
    anyons_star, anyons_plaq = lc.find_anyons_all()
# print("Anyons________>")
# print(anyons_star)
# print(anyon_plaq)
# print("-----------------<")

    sc.plot("star")

    match_star = matching.match(distance, anyons_star, topology,
                            "star", time=0, weights=weights)
    # match_plaq = matching.match(distance, anyons_plaq, topology,
    #                         "plaq", time=0, weights=weights)

# print("Matchings------->")
# print(match_star)

# print("Decoding now----->")
    sc.correct_error("star", match_star, time)
# sc.correct_error("plaq", match_plaq, time)


# sc.measure_stabilizer_type("star")
# sc.measure_stabilizer_type("plaq")

# if (sc.qubits[:, sc.tags == "Q"] == -1).any():
#     print("FAILURE CORRECTING")
# else:
#     print("SUCCESS CORRECTION")
    logical = sc.measure_logical()
    print(logical)

    sc.plot("star")

    plt.show()

if -1 in logical[0]:
    print("LOGICAL QUBIT ERR")

    # plt.close()
