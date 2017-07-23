"""
Surface code class.

created on: 19/07/17

@author: eduardo
"""

import numpy as np


class SurfaceCode:
    """
    Surface code class for either the toric or planar surfaces.

    STAR -- Q -- STAR -- Q -- STAR -- Q --
    |            |            |
    |            |            |
    Q    PLAQ    Q    PLAQ    Q    PLAQ
    |            |            |
    |            |            |
    PLAQ -- Q -- STAR -- Q -- STAR -- Q --
    |            |            |
    |            |            |
    Q    PLAQ    Q    PLAQ    Q    PLAQ
    |            |            |
    |            |            |
    PLAQ -- Q -- STAR -- Q -- STAR -- Q --
    |            |            |
    |            |            |

    Implements the operations:
        - Stabilizer measurement
        - Random noise
        - Noisy operations
        - Measurement errors

    All the information is saved on the (2, 2*size, 2*size) array qubits
    qubits[0] - X error and stabilizer measurements
    qubits[1] - Z error, stabilizer entries are not used
    """

    def __init__(self, size, surface="toroid"):
        """
        Init func.

        Parameters
        ----------
        size : int
            The base dimensions of the surface code.
        surface = "toroid" : string, optional
            Topology of the code.
        """
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
        self.qubits[1][self.tags != "Q"] = 9

        # transform to [[x1,x2,x3,...],[y1,y2,y3,...]]
        self.stars = self.stars.transpose()
        self.plaqs = self.plaqs.transpose()

        # Generate insterspersed stabilizer positions
        self.starsRound1 = self.stars[:, ::2]
        self.starsRound2 = self.stars[:, 1::2]
        self.plaqsRound1 = self.plaqs[:, ::2]
        self.plaqsRound2 = self.plaqs[:, 1::2]

    def measureStabilizer(self, pos, stabilizer, pNotComplete=0):
        """
        Measure stabilizers on the given position.

        Parameters
        ----------
        pos : array [[x1, x2, x3, ...], [y1, y2, y3, ...]]
            Positions of the stabilizers to be measured.
        stabilizer : string - "star" or "plaq"
            Kind of stabilizer to be measured.
        pNotComplete=0 : float
            Probaility to not complete stabilizer measurement.
        """
        _, c, t = self._selectStabilizer(stabilizer)

        if pNotComplete != 0:
            pos = self._incompleteMeasuerement(pos, pNotComplete)

        # TODO nearest indices func here: toroid planar
        t, b, l, r = self._stabilizerQubits(pos)
        self.qubits[0][pos[0], pos[1]] = (self.qubits[c][t[0], t[1]] *
                                          self.qubits[c][b[0], b[1]] *
                                          self.qubits[c][l[0], l[1]] *
                                          self.qubits[c][r[0], r[1]])

    def _incompleteMeasuerement(self, pos, pNotComplete):
        """Find stabilizers that are able to do a complete measurement."""
        # Calculate stabilizers that dont complete the measurement
        incomplete = (np.random.rand(len(pos)) < pNotComplete)
        # Remove them from the positions list
        newPos = np.delete(pos, np.where(incomplete), 1)
        return newPos

    def measureAllStabilizer(self, stabilizer):
        """
        Measure ALL stabilizers of a given type.

        Parameters
        ----------
        stabilizer: string - either "star" or "plaq"
        """
        pos, c, t = self._selectStabilizer(stabilizer)
        self.measureStabilizer(pos, stabilizer)

    def _stabilizerLie(self, tag, pLie):
        """Add measurement errot to a type of stabilizers."""
        # Add measurement error
        lie = 2*(np.random.rand(self.numberStabilizers) > pLie) - 1
        self.qubits[0][self.tags == tag] *= lie

    def _stabilizerLieAll(self, pLie):
        """Add measurement error to BOTH stars and plaqs stabilizers."""
        # Add measurement error to BOTH stabilizers
        lie = 2*(np.random.rand(2*self.numberStabilizers) > pLie) - 1
        self.qubits[0][self.tags != "Q"] *= lie

    def _stabilizerQubits(self, pos):
        """Find qubits corresponing to the given stabilizers."""
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

    def _applyNoiseQubit(self, pX, pZ):
        """Apply random error to the data qubits."""
        # Create the noise elements
        p = np.array([[pX], [pZ]])
        noise = 2*(np.random.rand(2, self.numberQubits) > p) - 1

        # Apply the noise
        self.qubits[:, self.tags == "Q"] *= noise

    def _applyOperationError(self, pos, error):
        """Apply operation error."""
        posQubit1, posQubit2 = self._twoRandStabQubits(pos)

        errMeasurement, errQubit1, errQubit2 = error
        # TODO care with errors shape
        # Errors have shape:
        # [[e1X, e2X, e3X, ...], [e1Z, e2Z, e3Z, ...]]
        self.qubits[:, posQubit1[0], posQubit1[1]] *= errQubit1
        self.qubits[:, posQubit2[0], posQubit2[1]] *= errQubit2

        # Apply error to stabilizer
        self.qubitts[0][pos[0], pos[1]] *= errMeasurement

    def NoisyMeasurement(self, stabilizer, errorVec, errorSum, pNotComplete=0):
        """
        Does a noisy measurement on the stabilizer type.
        The measurement is done in 2 rounds of interspersed stabilizers.
        """
        # Specify stabilizer
        pos1, pos2, c, t = self._selectStabilizerRounds(stabilizer)

        # Find all set of probabilistic errors
        err = np.random.rand(self.numberStabilizers)
        errIndex = np.zeros(self.numberStabilizers)
        for i in range(self.numberStabilizers):
            # Index of the actual error is the last True in less than rand number
            errIndex[i] = np.where(errorSum > err[i])[0][0]

        self.operationError(pos1, errorVec[errIndex[:len(pos1)]])
        self.measureStabilizer(pos1, stabilizer, pNotComplete)
        # TODO :or ?
        # self.measureStabilizer(pos, stabilizer)

        self.measureStabilizer(pos2, stabilizer, pNotComplete)
        self.operationError(pos2, errorVec[errIndex[len(pos1):]])


    def _twoRandStabQubits(self, pos):
        """Select to random corresponding qubits for each stabilizer."""
        options = np.array([[-1, 0], [1, 0], [0, -1], [0, 1]])
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
        """Return useful parameters for stabilizer type."""
        if stabilizer == "star":
            pos = self.stars
            c = 0
            t = "S"
        if stabilizer == "plaq":
            pos = self.plaqs
            c = 1
            t = "P"

        return pos, c, t

    def _selectStabilizerRounds(self, stabilizer):
        """Return useful parameters for stabilizer type."""
        if stabilizer == "star":
            pos1 = self.starsRound1
            pos2 = self.starsRound2
            c = 0
            t = "S"
        if stabilizer == "plaq":
            pos1 = self.plaqsRound1
            pos2 = self.plaqsRound2
            c = 1
            t = "P"

        return pos1, pos2, c, t

    def reset(self):
        sc.qubits.fill(1)
