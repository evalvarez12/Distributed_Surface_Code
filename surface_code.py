"""
Surface code class
bla bla
created on: 19/07/17

@author: eduardo
"""

import numpy as np


class SurfaceCode:
    """
    Surface code class
    """

    def __init__(self, size, surface="toroid"):
        self.size = size
        self.numberQubits = 2*size**2
        self.numberStabilizers = size**2
        self.surface = surface
        # Array with the quibits to mark erors
        # self.qubits[0] marks the X errors
        # self.qubits[1] marks the Z erros
        self.qubits = np.ones((2, 2*size, 2*size))

        # Positions of all stabilizers
        self.plaqs = []
        self.stars = []
        # Array with tags: Q, S or P useful for indexing
        self.tags = np.ones((2*size, 2*size), dtype=str)
        self.tags.fill("Q")
        for x in range(0, 2*size, 2):
            for y in range(0, 2*size, 2):
                self.tags[x, y] = "S"
                self.tags[x+1, y+1] = "P"

                self.stars += [[x, y]]
                self.plaqs += [[x+1, y+1]]
        self.stars = np.array(self.stars)
        self.plaqs = np.array(self.plaqs)

        # Fill the unused second entries of the stabilizers
        # with a 9 to mark
        self.qubits[1][ self.tags != "Q"] = 9

        # transform to [[x1,x2,x3,...],[y1,y2,y3,...]]
        self.stars = self.stars.transpose()
        self.plaqs = self.plaqs.transpose()

        # Generate insterspersed stabilizer positions
        self.starsRound1 = self.stars[:, ::2]
        self.starsRound2 = self.stars[:, 1::2]
        self.plaqsRound1 = self.plaqs[:, ::2]
        self.plaqsRound2 = self.plaqs[:, 1::2]

    def measureStabilizer(self, pos, stabilizer, pLie=0):
        """
        Measures all stabilizers
        """
        _, c, t = self._selectStabilizer(stabilizer)

        # TODO nearest indices func here: toroid planar
        t, b, l, r = self._stabilizerQubits(pos)
        self.qubits[0][pos[0], pos[1]] = (self.qubits[c][t[0], t[1]] *
                                          self.qubits[c][b[0], b[1]] *
                                          self.qubits[c][l[0], l[1]] *
                                          self.qubits[c][r[0], r[1]])



    def measureAllStabilizer(self, stabilizer, pLie=0):
        pos, c, t = self._selectStabilizer(stabilizer)
        self.measureStabilizer(pos, stabilizer, pLie)
        if pLie != 0:
            self.applyStabilizerLie(t, pLie)

    def applyStabilizerLie(self, tag, pLie):
        # Add measurement error
        lie = 2*(np.random.rand(self.numberStabilizers) > pLie) - 1
        self.qubits[0][self.tags == tag] *= lie

    def applyAllStabilizerLie(self, pLie):
        # Add measurement error to BOTH stabilizers
        lie = 2*(np.random.rand(2, self.numberStabilizers) > pLie) - 1
        self.qubits[0][self.tags != "Q"] *= lie




    def _stabilizerQubits(self, pos):
        top = pos + np.array([[-1], [0]])
        bottom = pos + np.array([[1], [0]])
        left = pos + np.array([[0], [-1]])
        right = pos + np.array([[0], [1]])

        if self.surface == "toroid":
            bottom = bottom % (2*self.size)
            right = right % (2*self.size)

        # TODO add this plane
        # if self.surface == "plane":

        return top, bottom, left, right


    def applyQubitErrors(self, pX, pZ):
        """
        Apply random error to the data qubits
        """
        # Create the noise elements
        p = np.array([[pX], [pZ]])
        noise = 2*(np.random.rand(2, self.numberQubits) > p) - 1

        # Apply the noise
        self.qubits[:, self.tags == "Q"] *= noise


    def applyNoisyMeasurement(self, stabilizer, errorVec, errorSum, notCompleteProb=0):
        # Specify stabilizer
        pos, c, t = self._selectStabilizer(stabilizer)

        # Find probabilistic error
        err = np.random.rand()
        # Index of the actual error is the last True in less than rand number
        errIndex = np.where(errorSum < err)[0][-1]


    # def applyErrors(self, pos, errors):


    def _twoStatilizerQubits(self, pos):
        options = np.array([[-1, 0],[1, 0],[0, -1], [0, 1]])
        # TODO optimize this for
        # Draw the selected options
        choices = np.array([np.random.choice([0, 1, 2, 3], 2, False) for i in range(len(pos[0]))])
        choices = choices.transpose()
        # print("-----------CHOICES--------")
        # print(choices)
        # print("--------------------------")
        displacement1 = options[choices[0]].transpose()
        displacement2 = options[choices[1]].transpose()

        # Position of the selected qubits
        a = pos + displacement1
        b = pos + displacement2

        # Adjust for surface type
        if self.surface == "toroid":
            a = a % (2*self.size)
            b = b % (2*self.size)

        # TODO add this plane
        # if self.surface == "plane":

        return a, b

    def _selectStabilizer(self, stabilizer):
        if stabilizer == "star":
            pos = self.stars
            c = 0
            t = "S"
        if stabilizer == "plaquette":
            pos = self.plaqs
            c = 1
            t = "P"

        return pos, c, t
