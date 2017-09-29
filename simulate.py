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


distance = 20
topology = "toroid"

sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)

sc._apply_noise_qubit(.1, .0)
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
print(sc.measure_logical())
sc.plot("star")
