"""
Layers class for 3d codes.

@author: eduardo
"""

import numpy as np

class Layers:
    """
    Layers to store the anyon positions in the surface codes.
    """

    def __init__(self, size):
        self.size = size
        self.numberStabilizers = size*size

        self.syndromeStar = np.ones(self.numberStabilizers)
        self.syndromePlaq = np.ones(self.numberStabilizers)

        self.parityStar = []
        self.parityPlaq = []

        self.stars = []
        self.plaqs = []
        for x in range(0, 2*size, 2):
            for y in range(0, 2*size, 2):
                self.stars += [[x, y]]
                self.plaqs += [[x+1, y+1]]

        self.stars = np.array(self.stars)
        self.plaqs = np.array(self.plaqs)

    def getTime(self):
        return len(self.parityStar)

    def add(self, codeStars, codePlaqs):
        # TODO redo this
        starLayer = codeStars
        plaqLayer = codePlaqs

        starLayer *= self.syndromeStar
        plaqLayer *= self.syndromePlaq

        self.syndromeStar = codeStars
        self.syndromePlaq = codePlaqs


        self.parityStar += [starLayer]
        self.parityPlaq += [plaqLayer]


    def findAnyons(self):
        # TODO improve all of this
        time = len(self.parityStar)

        anyonsStar = []
        anyonsPlaq = []
        for i in range(time):
            stars = self.stars[np.where(self.parityStar[i] == -1)]
            plaqs = self.plaqs[np.where(self.parityPlaq[i] == -1)]

            anyonsStar += [np.append(s, i) for s in stars]
            anyonsPlaq += [np.append(p, i) for p in plaqs]

        return np.array(anyonsStar), np.array(anyonsPlaq)
