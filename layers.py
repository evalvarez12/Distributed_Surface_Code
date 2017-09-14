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
        self.numberStabilizers = size*size

        # Lists to save the syndromes of previous layer
        self.pastSyndromeStar = np.ones(self.numberStabilizers)
        self.pastSyndromePlaq = np.ones(self.numberStabilizers)

        # Lists to save all the layers
        self.syndromesStar = []
        self.syndromesPlaq = []

        # The positions of the stars and plaqs for a code
        self.stars = []
        self.plaqs = []
        for x in range(0, 2*size, 2):
            for y in range(0, 2*size, 2):
                self.stars += [[x, y]]
                self.plaqs += [[x+1, y+1]]

        self.stars = np.array(self.stars)
        self.plaqs = np.array(self.plaqs)

    def getTime(self):
        return len(self.syndromesStar)

    def add(self, codeStars, codePlaqs):
        # New layer is obtained by comparing the previous one
        # with the new one, so physical errors are only carried once
        newStarLayer = codeStars * self.pastSyndromeStar
        newPlaqLayer = codePlaqs * self.pastSyndromePlaq

        # Save current code stabilizer measurements
        self.pastSyndromeStar = codeStars
        self.pastSyndromePlaq = codePlaqs

        # Add new layer of syndromes
        self.syndromesStar += [newStarLayer]
        self.syndromesPlaq += [newPlaqLayer]


    def findAnyons(self):
        # Get number of layers
        time = len(self.syndromesStar)

        anyonsStar = []
        anyonsPlaq = []
        # Go trough all layers finding the syndromes
        for i in range(time):
            stars = self.stars[np.where(self.syndromesStar[i] == -1)]
            plaqs = self.plaqs[np.where(self.syndromesPlaq[i] == -1)]

            # Save the physical and time positions of each syndrome
            anyonsStar += [np.append(s, i) for s in stars]
            anyonsPlaq += [np.append(p, i) for p in plaqs]

        return np.array(anyonsStar), np.array(anyonsPlaq)
