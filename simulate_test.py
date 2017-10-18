"""
File to test the surface code simulation.

created-on: 23/07/17
@author: eduardo
"""
import unittest
import numpy as np
import matplotlib.pyplot as plt
import surface_code
import layers
import matching

# Define the parameters
distance = 20
topology = "toric"
time_steps = 20
weights = [1, 1]

# Initialize objects
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)



class TestDecode(unittest.TestCase):

    def test_make_graph(self):
        for repetitions in range(100):
            # Perform measurements
            for i in range(time_steps):
                sc.apply_qubit_error(.1, .0)
                sc.measure_all_stablizers()
                sc._stabilizer_lie("S", .00)
                lc.add()


                time = lc.get_time()
                anyons_star, anyons_plaq = lc.find_anyons_all()

                match_star = matching.match(distance, anyons_star, topology,
                "star", time, weights=weights)
                match_plaq = matching.match(distance, anyons_plaq, topology,
                "plaq", time, weights=weights)

                sc.correct_error("star", match_star, time)
                sc.correct_error("plaq", match_plaq, time)

                logical = sc.measure_logical()
                sc.reset()
                lc.reset()
                print(logical)
                if -1 in logical[0]:
                    print("LOGICAL QUBIT ERR")


if __name__ == '__main__':
    unittest.main()
