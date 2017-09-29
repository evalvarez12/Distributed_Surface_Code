"""
Layers class for 3d codes.

@author: eduardo
"""

import numpy as np

class Layers:
    """
    Layers to store the anyon positions in the surface codes.
    Works as a 3d surface code where time is the third dimension.
    """

    def __init__(self, surface_code):
        # Save the suface code object
        self.surface_code = surface_code

        number_stabs = surface_code.number_stabs
        # Lists to save the syndromes of previous layer
        self.past_syndrome_star = np.ones(number_stabs)
        self.past_syndrome_plaq = np.ones(number_stabs)

        # Lists to save all the layers
        self.syndromes_star = []
        self.syndromes_plaq = []

        self.stars_pos = self.surface_code.stars.transpose()
        self.plaqs_pos = self.surface_code.plaqs.transpose()

    def get_time(self):
        return len(self.syndromes_star)

    def add(self):
        # New layer is obtained by comparing the previous one
        # with the new one, so physical errors are only carried once
        code_stars = self.surface_code.get_stars()
        code_plaqs = self.surface_code.get_plaqs()
        new_star_layer = code_stars * self.past_syndrome_star
        new_plaq_layer = code_plaqs * self.past_syndrome_plaq

        # Save current code stabilizer measurements
        self.past_syndrome_star = code_stars
        self.past_syndrome_plaq = code_plaqs

        # Add new layer of syndromes
        self.syndromes_star += [new_star_layer]
        self.syndromes_plaq += [new_plaq_layer]

    def find_anyons_all(self):
        anyons_star = self.find_anyons("star")
        anyons_plaq = self.find_anyons("plaq")
        return anyons_star, anyons_plaq

    def find_anyons(self, stabilizer):
        # Specify stabilizers
        if stabilizer == "plaq":
            stabs_pos = self.plaqs_pos
            syndromes = self.syndromes_plaq
        elif stabilizer == "star":
            stabs_pos = self.stars_pos
            syndromes = self.syndromes_star
        else:
            raise ValueError("Incorrect stabilizer argument")

        # Get number of layers
        time = len(syndromes)

        # Initialize anyons array, with flag empty
        anyons = np.array([0, 0, 0])
        empty = True

        # Go trough all layers finding the syndromes
        for t in range(time):
            indices = np.where(syndromes[t] == -1)[0]
            if len(indices) == 0:
                continue

            stabs = stabs_pos[indices]
            times = np.array([[t]] * len(indices))
            stabs = np.concatenate((stabs, times), 1)

            # Save the physical and time positions of each syndrome
            anyons = np.vstack((anyons, stabs))
            empty = False

        # Anyons array format is:
        # [[x1, y1, t1], [x2, y2, t2], [x3, y3, t3], ...]
        if empty:
            return []
        return anyons[1:]
