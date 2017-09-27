"""
File to test the surface code simulation.

created-on: 23/07/17
@author: eduardo
"""
import numpy as np

import surface_code
import layers
import matching


size = 3

sc = surface_code.SurfaceCode(size)
lc = layers.Layers(size)

sc._apply_noise_qubit(.2, 0)

sc.measure_stabilizer_type("star")
# sc.measure_stabilizer_type("plaq")

lc.add(sc.get_stars(), sc.get_plaqs())


anyon_star, anyon_plaq = lc.find_anyons_all()
print("Anyons________>")
print(anyon_star)
print("-----------------<")
# print(anyon_plaq)

match = matching.match_toric3D(size, anyon_star)
print("Matchings------->")
print(match)
print("Decoding now----->")
sc.correct_error("star", match)
