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
    bla bla
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

    def measure(self, stabilizer, pLie=0):
        """
        Measures all stabilizers
        """
        if stabilizer == "star":
            pos = self.stars
            c = 0
        if stabilizer == "plaquette":
            pos = self.plaqs
            c = 1

        # TODO nearest indices func here: toroid planar
        t, b, l, r = self._stabilizerQubits(pos)
        # print(pos)
        # print(t)
        # print(b)
        # print(l)
        # print(r)
        self.qubits[0][pos[0], pos[1]] = (self.qubits[c][t[0], t[1]] *
                                          self.qubits[c][b[0], b[1]] *
                                          self.qubits[c][l[0], l[1]] *
                                          self.qubits[c][r[0], r[1]])

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
        bla
        """
        # Create the noise elements
        p = np.array([[pX], [pZ]])

        noise = 2*(np.random.rand(2, self.numberQubits) > p) - 1

        # randX = 2*(np.random.rand(self.numberQubits) < pX) - 1
        # randZ = 2*(np.random.rand(self.numberQubits) < pZ) - 1
        # noise = np.array(randX, randZ)

        # Apply the noise
        self.qubits[:, self.tags == "Q"] *= noise
