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


distance = 10
topology = "toroid"
time_steps = 10
weights = [1, 2]
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)
for i in range(time_steps):
    sc.apply_qubit_error(.005, .0)
    sc.measure_all_stablizers()
    sc._stabilizer_lie("S", .009)
    lc.add()

# # sc.measure_all_stablizers()
# lc.add()
time = lc.get_time()

anyons_star, anyons_plaq = lc.find_anyons_all()
print("Anyons________>")
print(anyons_star)
# print(anyon_plaq)
print("-----------------<")
# print(anyon_plaq)

sc.plot("star")

if topology == "toroid":
    match_star = matching.match_toric_3D(distance, anyons_star, time, weights=weights)
    match_plaq = matching.match_toric_3D(distance, anyons_plaq, time,weights=weights)
else:
    match_star = matching.match_planar_3D(distance, anyons_star, "star", time, weights=weights)
    match_plaq = matching.match_planar_3D(distance, anyons_plaq, "plaq", time, weights=weights)


print("Matchings------->")
print(match_star)
print("Decoding now----->")
sc.correct_error("star", match_star, time)
sc.correct_error("plaq", match_plaq, time)


# sc.measure_stabilizer_type("star")
# sc.measure_stabilizer_type("plaq")

if (sc.qubits[:, sc.tags == "Q"] == -1).any():
    print("FAILURE CORRECTING")
else:
    print("SUCCESS CORRECTION")
logical = sc.measure_logical()
print(logical)

sc.plot("star")


if -1 in logical[0]:
    print("LOGICAL QUBIT ERR")

plt.show()
plt.close()
