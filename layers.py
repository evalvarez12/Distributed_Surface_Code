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

    def __init__(self, size):
        self.size = size
        self.number_stabilizers = size*size

        # Lists to save the syndromes of previous layer
        self.past_syndrome_star = np.ones(self.number_stabilizers)
        self.past_syndrome_plaq = np.ones(self.number_stabilizers)

        # Lists to save all the layers
        self.syndromes_star = []
        self.syndromes_plaq = []

        # The positions of the stars and plaqs for a code
        self.stars = []
        self.plaqs = []
        for x in range(0, 2*size, 2):
            for y in range(0, 2*size, 2):
                self.stars += [[x, y]]
                self.plaqs += [[x+1, y+1]]

        self.stars = np.array(self.stars)
        self.plaqs = np.array(self.plaqs)

    def get_time(self):
        return len(self.syndromes_star)

    def add(self, code_stars, code_plaqs):
        # New layer is obtained by comparing the previous one
        # with the new one, so physical errors are only carried once
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
            stabs = self.plaqs
            syndromes = self.syndromes_plaq
        elif stabilizer == "star":
            stabs = self.stars
            syndromes = self.syndromes_star
        else:
            raise ValueError("Incorrect stabilizer argument")

        # Get number of layers
        time = len(self.syndromes_star)

        anyons = []
        # Go trough all layers finding the syndromes
        for i in range(time):
            stabs = stabs[np.where(syndromes[i] == -1)]

            # Save the physical and time positions of each syndrome
            anyons += [np.append(s, i) for s in stabs]

        # Anyons array format is:
        # [[x1, y1, t1], [x2, y2, t2], [x3, y3, t3], ...]
        return np.array(anyons)
