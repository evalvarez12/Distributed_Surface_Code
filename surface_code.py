"""
Surface code class.

created on: 19/07/17

@author: eduardo
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


class SurfaceCode:
    """
    Surface code class for either the toric or planar surfaces.

    Q -- S -- Q -- S -- Q -- S --
         |         |         |
    P    Q    P    Q    P    Q
         |         |         |
    Q -- S -- Q -- S -- Q -- S --
         |         |         |
    P    Q    P    Q    P    Q
         |         |         |
    Q -- S -- Q -- S -- Q -- S --
         |         |         |
    P    Q    P    Q    P    Q
         |         |         |

    Implements the operations:
        - Stabilizer measurement
        - Random noise
        - Noisy operations
        - Measurement errors

    All the information is saved on the (2, 2*distance, 2*distance) array qubits
    qubits[0] - X error and stabilizer measurements
    qubits[1] - Z error, stabilizer entries are not used
    """

    def __init__(self, distance, surface):
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
        if self.surface == "toroid":
            self.number_data_qubits = 2*distance**2
            self.number_stabs = distance**2
            self.side = 2*distance

        elif self.surface == "planar":
            self.number_data_qubits = distance**2 + (distance - 1)**2
            self.number_stabs = (distance - 1)*distance
            self.side = 2*distance - 1

        ind1 = np.arange(1, self.side, 2)
        ind2 = np.arange(0, self.side, 2)

        # Array with the quibits to mark erors
        # self.qubits[0] marks the X errors
        # self.qubits[1] marks the Z erros

        starsy, starsx = np.meshgrid(ind1, ind2)
        plaqsy, plaqsx = np.meshgrid(ind2, ind1)

        self.qubits = np.ones((2, self.side, self.side))

        self.stars = np.vstack((starsx.flatten(), starsy.flatten()))
        self.plaqs = np.vstack((plaqsx.flatten(), plaqsy.flatten()))

        # Array with tags: Q, S or P useful for indexing
        self.tags = np.ones((self.side, self.side), dtype=str)
        self.tags.fill("Q")
        self.tags[self.stars[0], self.stars[1]] = "S"
        self.tags[self.plaqs[0], self.plaqs[1]] = "P"

        # Fill the unused second entries of the stabilizers
        # with a 9 to mark
        self.qubits[1][self.tags != "Q"] = 0

        # Generate insterspersed stabilizer positions
        self.stars_round1 = self.stars[:, ::2]
        self.stars_round2 = self.stars[:, 1::2]
        self.plaqs_round1 = self.plaqs[:, ::2]
        self.plaqs_round2 = self.plaqs[:, 1::2]

        # Color map stuff for plot
        self.cmap = colors.ListedColormap(['red', 'orange', 'white', 'green'])
        bounds = [-2.5, -1.5, 0, 1.5, 2.5]
        self.cmap_norm = colors.BoundaryNorm(bounds, self.cmap.N)

        # NOTE DEPRECATED
        # for x in range(0, 2*distance, 2):
        #     for y in range(0, 2*distance, 2):
        #         self.tags[x, y] = "S"
        #         self.tags[x+1, y+1] = "P"
        #
        #         self.stars += [[x, y]]
        #         self.plaqs += [[x+1, y+1]]
        # self.stars = np.array(self.stars)
        # self.plaqs = np.array(self.plaqs)


        # # transform to [[x1,x2,x3,...],[y1,y2,y3,...]]
        # self.stars = self.stars.transpose()
        # self.plaqs = self.plaqs.transpose()


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
        # TODO replace this with fancy indexing
        t, b, l, r = self._stabilizer_qubits(pos)
        vals_t = self.qubits[c][t[0], t[1]]
        vals_b = self.qubits[c][b[0], b[1]]
        vals_l = self.qubits[c][l[0], l[1]]
        vals_r = self.qubits[c][r[0], r[1]]

        if self.surface == "planar":
            # Invalidate the over border "results"
            # NOTE this needs to change
            vals_t[np.where(t[0] == -1)] = 1
            vals_b[np.where(b[0] == -1)] = 1
            vals_l[np.where(l[1] == -1)] = 1
            vals_r[np.where(r[1] == -1)] = 1

        self.qubits[0][pos[0], pos[1]] = (vals_t
                                          * vals_b
                                          * vals_l
                                          * vals_r)
        # NOTE DEPRECATED!
        # self.qubits[0][pos[0], pos[1]] = (self.qubits[c][t[0], t[1]] *
        #                                   self.qubits[c][b[0], b[1]] *
        #                                   self.qubits[c][l[0], l[1]] *
        #                                   self.qubits[c][r[0], r[1]])

    def _stabilizer_qubits(self, pos):
        """Find qubits corresponing to the given stabilizers."""
        top = pos + np.array([[-1], [0]])
        bottom = pos + np.array([[1], [0]])
        left = pos + np.array([[0], [-1]])
        right = pos + np.array([[0], [1]])

        if self.surface == "toroid":
            # Take the mod to account for cyclic boundaries
            bottom = bottom % self.side
            right = right % self.side
            # Top and left are automatically accounted for
            # when -1 is the index
        elif self.surface == "planar":
            # Mark the list indices which go outside of the lattice
            # with -1 for latter removal
            bottom[bottom >= self.side] = -1
            right[right >= self.side] = -1
            top[top <= -1] = -1
            left[left <= -1] = -1

        return top, bottom, left, right

    def measure_stabilizer_type(self, stabilizer, p_not_complete=0):
        """
        Measure ALL stabilizers of a given type.

        Parameters
        ----------
        stabilizer: string - either "star" or "plaq"
        """
        pos, c, t = self._select_stabilizer(stabilizer)
        self.measure_stabilizer(pos, c, p_not_complete)

    def _incomplete_measuerement(self, pos, p_not_complete):
        """Find stabilizers that are able to do a complete measurement."""
        # Calculate stabilizers that dont complete the measurement
        incomplete = (np.random.rand(len(pos)) < p_not_complete)
        # Remove them from the positions list
        new_pos = np.delete(pos, np.where(incomplete), 1)
        return new_pos

    def _stabilizer_lie(self, tag, p_lie):
        """Add measurement errot to a type of stabilizers."""
        # Add measurement error
        lie = 2*(np.random.rand(self.number_stabs) > p_lie) - 1
        self.qubits[0][self.tags == tag] *= lie

    def _stabilizer_lie_all(self, p_lie):
        """Add measurement error to BOTH stars and plaqs stabilizers."""
        # Add measurement error to BOTH stabilizers
        self._stabilizer_lie("S", p_lie)
        self._stabilizer_lie("P", p_lie)

    def _apply_noise_qubit(self, pX, pZ):
        """Apply random error to the data qubits."""
        # Create the noise elements
        p = np.array([[pX], [pZ]])
        noise = 2*(np.random.rand(2, self.number_data_qubits) > p) - 1

        # Apply the noise
        self.qubits[:, self.tags == "Q"] *= noise

    def _apply_operation_error(self, pos, error):
        """Apply operation error."""
        # NOTE This will change
        pos_qubit1, pos_qubit2 = self.two_rand_stab_qubits(pos)

        err_measurement, err_qubit1, err_qubit2 = error
        # TODO care with errors shape
        # Errors have shape:
        # [[e1X, e2X, e3X, ...], [e1Z, e2Z, e3Z, ...]]
        self.qubits[:, pos_qubit1[0], pos_qubit1[1]] *= err_qubit1
        self.qubits[:, pos_qubit2[0], pos_qubit2[1]] *= err_qubit2

        # Apply error to stabilizer
        self.qubitts[0][pos[0], pos[1]] *= err_measurement

    def operation_error(self, pos, stabilizer):
        pos, measurement_err, qubit_err = self.errors.get_errors(pos,
                                                                 stabilizer)

        # NOTE Working here
        self.qubits[:, pos_qubit1[0], pos_qubit1[1]] *= err_qubit1
        self.qubits[:, pos_qubit2[0], pos_qubit2[1]] *= err_qubit2


        # Apply error to stabilizer
        self.qubitts[0][pos[0], pos[1]] *= err_measurement



    def noisy_measurement(self, stabilizer, error_vec, error_sum, p_not_complete=0):
        """
        Does a noisy measurement on the stabilizer type.
        The measurement is done in 2 rounds of interspersed stabilizers.
        """
        # NOTE this will change
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


    # NOTE: This will be DEPRECATED
    # def _two_rand_stab_qubits(self, pos):
    #     """Select to random corresponding qubits for each stabilizer."""
    #     options = np.array([[-1, 0], [1, 0], [0, -1], [0, 1]])
    #     # TODO optimize this for
    #     # Draw the selected options
    #     choices = np.array([np.random.choice([0, 1, 2, 3], 2, False) for i in range(len(pos[0]))])
    #     choices = choices.transpose()
    #     # print("-----------CHOICES--------")
    #     # print(choices)
    #     # print("--------------------------")
    #     displacement1 = options[choices[0]].transpose()
    #     displacement2 = options[choices[1]].transpose()
    #
    #     # Position of the selected qubits
    #     a = pos + displacement1
    #     b = pos + displacement2
    #
    #     # Adjust for surface type
    #     if self.surface == "toroid":
    #         a = a % (2*self.distance)
    #         b = b % (2*self.distance)
    #
    #     # TODO add this plane
    #     # if self.surface == "plane":
    #
    #     return a, b

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
        """Parameter for stabilizer type when using interspersed rounds."""
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

    def plot(self, stabilizer):
        """Plot the surface code."""
        if stabilizer == "star":
            data = self.qubits[0].copy()
            data[self.tags == "S"] *= 2
            data[self.tags == "P"] = 1
        if stabilizer == "plaq":
            data = self.qubits[0].copy()
            data[self.tags == "Q"] = self.qubits[1][self.tags == "Q"]
            data[self.tags == "P"] *= 2
            data[self.tags == "S"] = 1

        plt.imshow(data, cmap=self.cmap, norm=self.cmap_norm)
        # plt.colorbar()
        plt.show()

    def reset(self):
        """Reset surface code to default configuration."""
        self.qubits.fill(1)
        self.qubits[1][self.tags != "Q"] = 0

    def get_stars(self):
        """Get the values of  all star qubits."""
        return self.qubits[0][self.tags == "S"]

    def get_plaqs(self):
        """Get the values of all plaq qubits."""
        return self.qubits[0][self.tags == "P"]

    def measure_logical(self):
        """Meausure the logical qubits."""
        if self.surface == "toroid":
            # TODO check here
            # X1 second row - Z1 first row
            X1 = np.prod(self.qubits[0, 1, 1::2])
            Z1 = np.prod(self.qubits[1, 0, 0::2])

            # X2 first column - Z1 second column
            X2 = np.prod(self.qubits[0, 0::2, 0])
            Z2 = np.prod(self.qubits[1, 1::2, 1])
            X = [X1, X2]
            Z = [Z1, Z2]
        elif self.surface == "planar":
            # X first column - Z first row
            X = [np.prod(self.qubits[0, 0::2, 0])]
            Z = [np.prod(self.qubits[1, 0, 0::2])]
        return X, Z


    def correct_error(self, error_type, match):
        _, c, t = self._select_stabilizer(error_type)

        #Taking the transpose of match gives:
        # [[[p1x, p2x, p3x, ...], [q1x, q2x, q3x, ...]],
        #  [[p1y, p2y, p3y, ...], [q1y, q2y, q3y, ...]],
        #  [[p1t, p2t, p3t, ...], [q1t, q2t, q3t, ...]]]
        match = match.transpose()

        # Distances in each coordinate
        dx = np.abs(match[0][0] - match[0][1])
        dy = np.abs(match[1][0] - match[1][1])
        dt = np.abs(match[2][0] - match[2][1])


        if self.surface == "planar":
            for i in range(len(dx)):
                # Start and end points of the path
                startx = match[0][0][i]
                endx = match[0][1][i]
                starty = match[1][0][i]
                # Create steps and join them
                stepsx = np.arange(1, dx[i], 2) + startx
                stepsy = np.arange(1, dy[i], 2) + starty
                stepsx = np.append(stepsx, [endx] * len(stepsy))
                stepsy = np.append([starty] * len(stepsx), stepsy)

                # Apply error correction path
                self.qubits[c][stepsx, stepsy] *= -1

        elif self.surface == "toroid":
            for i in range(len(dx)):

                if dx[i] < self.distance:
                    # "Normal" path in the surface
                    startx = match[0][0][i]
                    endx = match[0][1][i]
                else:
                    # Path in the far side of the torus
                    d = dx[i] % self.distance
                    startx = match[0][1][i]
                    endx = match[0][0][i]

                stepsx = np.arange(1, dx[i], 2) + startx
                stepsx = stepsx % self.side

                if dy[i] < self.distance:
                    starty = match[1][0][i]
                else:
                    d = dy[i] % self.distance
                    starty = match[1][1][i]
                stepsy = np.arange(1, dy[i], 2) + starty
                stepsy = stepsy % self.side
                # Join all stepss
                stepsx = np.append(stepsx, [endx] * len(stepsy))
                stepsy = np.append([starty] * len(stepsx), stepsy)

                # Apply error correction path
                self.qubits[c][stepsx, stepsy] *= -1
