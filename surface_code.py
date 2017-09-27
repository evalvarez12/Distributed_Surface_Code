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
    STAR -- Q -- STAR -- Q -- STAR -- Q --
    |            |            |
    |            |            |
    Q    PLAQ    Q    PLAQ    Q    PLAQ
    |            |            |
    |            |            |
    STAR -- Q -- STAR -- Q -- STAR -- Q --
    |            |            |
    |            |            |

    Implements the operations:
        - Stabilizer measurement
        - Random noise
        - Noisy operations
        - Measurement errors

    All the information is saved on the (2, 2*distance, 2*distance) array qubits
    qubits[0] - X error and stabilizer measurements
    qubits[1] - Z error, stabilizer entries are not used
    """

    def __init__(self, distance, surface="toroid"):
        """
        Init func.

        Parameters
        ----------
        distance : int
            The distance of the surface of the surface code.
        surface = "toroid" : string, optional
            Topology of the code.
        """
        self.distance = distance
        self.surface = surface
        # TODO: this is for toroid check for plannar
        self.number_qubits = 2*distance**2
        self.number_stabilizers = distance**2
        # Array with the quibits to mark erors
        # self.qubits[0] marks the X errors
        # self.qubits[1] marks the Z erros
        self.qubits = np.ones((2, 2*distance, 2*distance))

        # Positions of all stabilizers
        self.plaqs = []
        self.stars = []
        # Array with tags: Q, S or P useful for indexing
        self.tags = np.ones((2*distance, 2*distance), dtype=str)
        self.tags.fill("Q")
        for x in range(0, 2*distance, 2):
            for y in range(0, 2*distance, 2):
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
        self.stars_round1 = self.stars[:, ::2]
        self.stars_round2 = self.stars[:, 1::2]
        self.plaqs_round1 = self.plaqs[:, ::2]
        self.plaqs_round2 = self.plaqs[:, 1::2]

    def measure_stabilizer(self, pos, c, p_not_complete=0):
        """
        Measure stabilizers on the given position.

        Parameters
        ----------
        pos : array [[x1, x2, x3, ...], [y1, y2, y3, ...]]
            Positions of the stabilizers to be measured.
        c   : 0 or 1 channel of the stabilizer type
        stabilizer : string - "star" or "plaq"
            Kind of stabilizer to be measured.
        p_not_complete=0 : float
            Probaility to not complete stabilizer measurement.
        """

        if p_not_complete != 0:
            pos = self._incomplete_measuerement(pos, p_not_complete)

        # TODO nearest indices func here: toroid planar
        t, b, l, r = self._stabilizer_qubits(pos)
        self.qubits[0][pos[0], pos[1]] = (self.qubits[c][t[0], t[1]] *
                                          self.qubits[c][b[0], b[1]] *
                                          self.qubits[c][l[0], l[1]] *
                                          self.qubits[c][r[0], r[1]])

    def _incomplete_measuerement(self, pos, p_not_complete):
        """Find stabilizers that are able to do a complete measurement."""
        # Calculate stabilizers that dont complete the measurement
        incomplete = (np.random.rand(len(pos)) < p_not_complete)
        # Remove them from the positions list
        new_pos = np.delete(pos, np.where(incomplete), 1)
        return new_pos

    def measure_stabilizer_type(self, stabilizer, p_not_complete=0):
        """
        Measure ALL stabilizers of a given type.

        Parameters
        ----------
        stabilizer: string - either "star" or "plaq"
        """
        pos, c, t = self._select_stabilizer(stabilizer)
        self.measure_stabilizer(pos, c, p_not_complete)

    def _stabilizer_lie(self, tag, p_lie):
        """Add measurement errot to a type of stabilizers."""
        # Add measurement error
        lie = 2*(np.random.rand(self.number_stabilizers) > p_lie) - 1
        self.qubits[0][self.tags == tag] *= lie

    def _stabilizer_lie_all(self, p_lie):
        """Add measurement error to BOTH stars and plaqs stabilizers."""
        # Add measurement error to BOTH stabilizers
        lie = 2*(np.random.rand(2*self.number_stabilizers) > p_lie) - 1
        self.qubits[0][self.tags != "Q"] *= lie

    def _stabilizer_qubits(self, pos):
        """Find qubits corresponing to the given stabilizers."""
        top = pos + np.array([[-1], [0]])
        bottom = pos + np.array([[1], [0]])
        left = pos + np.array([[0], [-1]])
        right = pos + np.array([[0], [1]])

        if self.surface == "toroid":
            bottom = bottom % (2*self.distance)
            right = right % (2*self.distance)

        # TODO add this plane
        # if self.surface == "plane":

        return top, bottom, left, right

    def _apply_noise_qubit(self, pX, pZ):
        """Apply random error to the data qubits."""
        # Create the noise elements
        p = np.array([[pX], [pZ]])
        noise = 2*(np.random.rand(2, self.number_qubits) > p) - 1

        # Apply the noise
        self.qubits[:, self.tags == "Q"] *= noise

    def _apply_operation_error(self, pos, error):
        """Apply operation error."""
        pos_qubit1, pos_qubit2 = self.two_rand_stab_qubits(pos)

        err_measurement, err_qubit1, err_qubit2 = error
        # TODO care with errors shape
        # Errors have shape:
        # [[e1X, e2X, e3X, ...], [e1Z, e2Z, e3Z, ...]]
        self.qubits[:, pos_qubit1[0], pos_qubit1[1]] *= err_qubit1
        self.qubits[:, pos_qubit2[0], pos_qubit2[1]] *= err_qubit2

        # Apply error to stabilizer
        self.qubitts[0][pos[0], pos[1]] *= err_measurement

    def noisy_measurement(self, stabilizer, error_vec, error_sum, p_not_complete=0):
        """
        Does a noisy measurement on the stabilizer type.
        The measurement is done in 2 rounds of interspersed stabilizers.
        """
        # Specify stabilizer
        pos1, pos2, c, t = self._select_stabilizerRounds(stabilizer)

        # Find all set of probabilistic errors
        err = np.random.rand(self.number_stabilizers)
        err_index = np.zeros(self.number_stabilizers)
        for i in range(self.number_stabilizers):
            # Index of the actual error is the last True in less than rand number
            err_index[i] = np.where(error_sum > err[i])[0][0]

        self.operation_error(pos1, error_vec[err_index[:len(pos1)]])
        self.measure_stabilizer(pos1, c, p_not_complete)
        # TODO :or ?
        # self.measureStabilizer(pos, c)

        self.measure_stabilizer(pos2, c, p_not_complete)
        self.operation_error(pos2, error_vec[err_index[len(pos1):]])


    def _two_rand_stab_qubits(self, pos):
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
            a = a % (2*self.distance)
            b = b % (2*self.distance)

        # TODO add this plane
        # if self.surface == "plane":

        return a, b

    def _select_stabilizer(self, stabilizer):
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

    def _select_stabilizer_rounds(self, stabilizer):
        """Parameters for stabilizer type when using interspersed rounds."""
        if stabilizer == "star":
            pos1 = self.stars_round1
            pos2 = self.stars_round2
            c = 0
            t = "S"
        if stabilizer == "plaq":
            pos1 = self.plaqs_round1
            pos2 = self.plaqs_round2
            c = 1
            t = "P"

        return pos1, pos2, c, t

    def reset(self):
        """Reset surface code to default configuration."""
        self.qubits.fill(1)
        self.qubits[1][self.tags != "Q"] = 9

    def get_stars(self):
        """Get the position of  all star qubits."""
        return self.qubits[0][self.tags == "S"]

    def get_plaqs(self):
        """Get the positiono of all plaq qubits."""
        return self.qubits[0][self.tags == "P"]

    def measure_logical(self):
        """Meausure the logical qubits."""
        # X1 first row - Z1 second row
        X1 = np.prod(self.qubits[0, 0, 1::2])
        Z1 = np.prod(self.qubits[1, 1, 0::2])

        # X2 first column - Z1 second column
        X2 = np.prod(self.qubits[0, 1::2, 0])
        Z2 = np.prod(self.qubits[1, 0::2, 1])

        return [[X1, X2], [Z1, Z2]]


    def correct_error(self, error_type, match):
        # TODO redo all  of this?
        _, c, t = self._select_stabilizer(error_type)


        m = 2*self.distance
        for pair in match:
            print("Pair:", pair)
            px, py, _ = pair[0]
            qx, qy, _ = pair[1]

            dx = (qx - px) % m
            dy = (qy - py) % m

            if dx < m - dx:
                endx = qx
                stepsx = np.arange(1, dx, 2)
                stepsx = (stepsx + px) % m
                coord_y = np.ones_like(stepsx)*py
            else:
                endx = px
                stepsx = np.arange(1, m - dx, 2)
                stepsx = (stepsx + qx) % m
                coord_y = np.ones_like(stepsx)*qy

            if dy < m - dy:
                stepsy = (np.arange(1, dy, 2) + py) % m

            else:
                stepsy = (np.arange(1, m - dy, 2) + qy) % m

            stepsx = np.append(stepsx, [endx] * len(stepsy)).astype(int)
            stepsy = np.append(coord_y, stepsy).astype(int)

            print(stepsx)
            print(stepsy)
            self.qubits[c][stepsx, stepsy] *= -1
