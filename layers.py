"""
Layers class for holding the measurement sheets when performing imperfect measurements
on the surface code.

author: Eduardo Villasenor
created-on: 12/07/17
"""

import numpy as np
import matching

class Layers:
    """
    Layers to store the syndrome measurements in the surface codes.
    Works as a 3d surface code where time is the third dimension.
    """

    def __init__(self, surface_code):
        """
        Init function

        Parameters
        -----------
        surface_code : (SurfaceCode) surface code object
        """
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
        """Returns the number of measurements done so far"""
        return len(self.syndromes_star)

    def reset(self):
        """Erases all saved syndroms"""
        self.past_syndrome_star = np.ones(self.surface_code.number_stabs)
        self.past_syndrome_plaq = np.ones(self.surface_code.number_stabs)

        self.syndromes_star = []
        self.syndromes_plaq = []

    def decode(self, weights=[1, 1]):
        """Decode the measured syndromes using MWPM"""
        # Get anyons
        anyons_star, anyons_plaq = self.find_anyons_all()

        # The number of measurements done so far
        time = len(self.syndromes_star)
        # Decode the stabilizer measuements
        match_star = matching.match_simple(self.surface_code.distance,
                                           anyons_star,
                                           self.surface_code.surface,
                                           "star", weights=weights)
        match_plaq = matching.match_simple(self.surface_code.distance,
                                           anyons_plaq,
                                           self.surface_code.surface,
                                           "plaq", weights=weights)
        # Apply corrections to the surface code
        self.surface_code.correct_error("star", match_star)
        self.surface_code.correct_error("plaq", match_plaq)


    def add(self):
        """Add the current measuement status as another layer."""
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
        """Return all the stabilizer measuements so far"""
        anyons_star = self.find_anyons("star")
        anyons_plaq = self.find_anyons("plaq")
        return anyons_star, anyons_plaq

    def find_anyons(self, stabilizer):
        # Return a list with all the stabilizer measurements
        # Anyons array format is:
        # [[x1, y1, t1], [x2, y2, t2], [x3, y3, t3], ...]
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

        if empty:
            return []
        return anyons[1:]
