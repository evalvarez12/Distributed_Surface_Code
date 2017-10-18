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
distance = 20
topology = "toric"
time_steps = 1
weights = [1, 1]

# Initialize objects
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)

# Perform measurements
for i in range(time_steps):
    sc.apply_qubit_error(.1, .0)
    sc.measure_all_stablizers()
    sc._stabilizer_lie("S", .00)
    lc.add()


time = lc.get_time()
anyons_star, anyons_plaq = lc.find_anyons_all()
print("Anyons________>")
print(anyons_star)
# print(anyon_plaq)
print("-----------------<")

sc.plot("star")

match_star = matching.match(distance, anyons_star, topology,
                            "star", time=0, weights=weights)
match_plaq = matching.match(distance, anyons_plaq, topology,
                            "plaq", time=0, weights=weights)

print("Matchings------->")
print(match_star)

print("Decoding now----->")
sc.correct_error("star", match_star, time)
sc.correct_error("plaq", match_plaq, time)


# sc.measure_stabilizer_type("star")
# sc.measure_stabilizer_type("plaq")

# if (sc.qubits[:, sc.tags == "Q"] == -1).any():
#     print("FAILURE CORRECTING")
# else:
#     print("SUCCESS CORRECTION")
logical = sc.measure_logical()
print(logical)

sc.plot("star")


if -1 in logical[0]:
    print("LOGICAL QUBIT ERR")

plt.show()
plt.close()
