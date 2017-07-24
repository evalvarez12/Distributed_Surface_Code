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

sc._applyNoiseQubit(.5,.5)

sc.measureAllStabilizer("star")
sc.measureAllStabilizer("plaq")

lc.add(sc.getStars(), sc.getPlaqs())


anyonStar, anyonPlaq = lc.findAnyons()
print(anyonStar)
print(anyonPlaq)

match = matching.matchToric3D(size, anyonStar)
print(match)
