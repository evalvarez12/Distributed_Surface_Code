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


distance = 5
topology = "planar"

sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)
for i in range(1):
    sc._apply_noise_qubit(.09, .0)
    sc.measure_stabilizer_type("star")
    sc.measure_stabilizer_type("plaq")

    sc.plot("star")

    lc.add()

    anyon_star, anyon_plaq = lc.find_anyons_all()
    print("Anyons________>")
    print(anyon_star)
    print(anyon_plaq)
    print("-----------------<")
    # print(anyon_plaq)

    match_star = matching.match_toric3D(distance, anyon_star)
    match_plaq = matching.match_toric3D(distance, anyon_plaq)

    print("Matchings------->")
    # print(match)
    print("Decoding now----->")
    sc.correct_error("star", match_star)
    sc.correct_error("plaq", match_plaq)
    sc.measure_stabilizer_type("star")
    sc.measure_stabilizer_type("plaq")

    if (sc.qubits[0][sc.tags != "Q"] == -1).any():
        print("FAILURE CORRECTING")
    else:
        print("SUCCESS CORRECTION")
    logical = sc.measure_logical()
    print(logical)
    sc.plot("star")
    if -1 in logical[0]:
        print("LOGICAL QUBIT ERR")
        break
