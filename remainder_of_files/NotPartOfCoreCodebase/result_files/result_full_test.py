"""
Simple simulation to test the surface code simulation is working
without looking at plots.

created-on: 09/12/17
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt
import surface_code
import layers
import matching

"""
Test the planar code.
Perfect syndrome extraction.
"""

# Define the parameters
distance = 10
topology = "planar"
weights = [1, 1]
iterations = 5
p = 0.05
fail_rate_planar = 0

# Initialize objects
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)

# Perform measurements
for i in range(iterations):

    # Errors and measurements
    sc.apply_qubit_error(p, p)
    sc.measure_all_stabilizers()
    lc.add()

    # Get anyons
    anyons_star, anyons_plaq = lc.find_anyons_all()

    # Decode
    match_star = matching.match(distance, anyons_star, topology,
                                "star", time=0, weights=weights)
    match_plaq = matching.match(distance, anyons_plaq, topology,
                                "plaq", time=0, weights=weights)

    # Apply corrections
    sc.correct_error("star", match_star)
    sc.correct_error("plaq", match_plaq)

    # Check for errors in decoding and correcting
    sc.measure_stabilizer_type("star")
    sc.measure_stabilizer_type("plaq")
    if (sc.qubits[:, sc.tags != "Q"] == -1).any():
        print("FAILURE CORRECTING")

    logical = sc.measure_logical()

    if -1 in logical[0] or -1 in logical[1]:
        fail_rate_planar += 1

    lc.reset()
    sc.reset()


"""
Test the toric code.
Imperfect syndrome extraction.
"""

# Define the parameters
distance = 10
topology = "toric"
weights = [1, 1]
iterations = 5
p = 0.005
q = 0.005
cycles = 50
fail_rate_toric = 0

# Initialize objects
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)

# Perform measurements
for i in range(iterations):

    # Errors and measurements
    for s in range(cycles):
        sc.apply_qubit_error(p, p)
        sc.measure_all_stablizers()
        sc.apply_measurement_error(q)
        lc.add()

    # Get anyons
    anyons_star, anyons_plaq = lc.find_anyons_all()

    # Decode
    match_star = matching.match(distance, anyons_star, topology,
                                "star", time=cycles, weights=weights)
    match_plaq = matching.match(distance, anyons_plaq, topology,
                                "plaq", time=cycles, weights=weights)

    # Apply corrections
    sc.correct_error("star", match_star, cycles)
    sc.correct_error("plaq", match_plaq, cycles)

    # Round of perfect detection to eliminate stray errors
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

    logical = sc.measure_logical()

    if -1 in logical[0] or -1 in logical[1]:
        fail_rate_toric += 1

    lc.reset()
    sc.reset()

fail_rate_toric = fail_rate_toric/float(iterations)
fail_rate_planar = fail_rate_planar/float(iterations)
print("SUCCESS RATE TORIC: ", 1 - fail_rate_toric)
print("SUCCESS RATE PLANAR: ", 1 - fail_rate_planar)
