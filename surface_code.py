"""
Surface code class.

created on: 19/07/17

@author: eduardo
"""

import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import colors
import errors


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
        surface = "toric" / "planar" : string, optional
            Topology of the code.
        """
        if surface != "toric" and surface != "planar":
            raise ValueError("Incorrect surface argument SurfaceCode")
        self.distance = distance
        self.surface = surface
        # TODO: this is for toric check for plannar
        if self.surface == "toric":
            self.number_data_qubits = 2*distance**2
            self.number_stabs = distance**2
            self.side = 2*distance

        elif self.surface == "planar":
            self.number_data_qubits = distance**2 + (distance - 1)**2
            self.number_stabs = (distance - 1)*distance
            self.side = 2*distance - 1

        ind1 = np.arange(1, self.side, 2)
        ind2 = np.arange(0, self.side, 2)

        # Array with the qubits to mark erors
        # self.qubits[0] marks the X errors
        # self.qubits[1] marks the Z erros
        self.qubits = np.ones((2, self.side, self.side))

        starsy, starsx = np.meshgrid(ind1, ind2)
        plaqsy, plaqsx = np.meshgrid(ind2, ind1)
        self.stars = np.vstack((starsx.flatten(), starsy.flatten()))
        self.plaqs = np.vstack((plaqsx.flatten(), plaqsy.flatten()))

        # Array with tags: Q, S or P useful for indexing
        self.tags = np.ones((self.side, self.side), dtype=str)
        self.tags.fill("Q")
        self.tags[self.stars[0], self.stars[1]] = "S"
        self.tags[self.plaqs[0], self.plaqs[1]] = "P"
        # Fill the unused second entries of the stabilizers
        # with a 9 to mark
        self.qubits[1, self.tags != "Q"] = 0

        # Plane tags to mark boundary stabilizers
        # Only for planar topology
        if self.surface == "planar":
            self.plane = np.ones((self.side, self.side), dtype=str)
            self.plane.fill("o")
            self.plane[0] = "t"
            self.plane[-1] = "b"
            self.plane[:, 0] = "l"
            self.plane[:, -1] = "r"
            self.plane[self.tags == "Q"] = "-"

        # Generate insterspersed stabilizer positions
        self.stars_round1 = self.stars[:, ::2]
        self.stars_round2 = self.stars[:, 1::2]
        self.plaqs_round1 = self.plaqs[:, ::2]
        self.plaqs_round2 = self.plaqs[:, 1::2]

        # Color map stuff for plot
        # self.cmap = colors.ListedColormap(['red', 'orange', 'white', 'green'])
        # bounds = [-2.5, -1.5, 0, 1.5, 2.5]
        # self.cmap_norm = colors.BoundaryNorm(bounds, self.cmap.N)


    def init_error_obj(self, surface, ps, pm, pg, eta, a0, a1, theta, protocol):
        self.errors = errors.Generator(surface=self.surface, ps=ps, pm=pm,
                                       pg=pg, eta=eta, a0=a0, a1=a1,
                                       theta=theta, protocol=protocol)

    def measure_all_stablizers(self, p_not_complete=0):
        """Measure all stabilizer in the code."""
        self.measure_stabilizer_type("star", p_not_complete)
        self.measure_stabilizer_type("plaq", p_not_complete)

    def measure_stabilizer_type(self, stabilizer, p_not_complete=0):
        """
        Measure ALL stabilizers of a given type.

        Parameters
        ----------
        stabilizer: string - either "star" or "plaq"
        """
        if stabilizer == "star":
            pos = self.stars
            c = 0
        if stabilizer == "plaq":
            pos = self.plaqs
            c = 1

        self.measure_stabilizer(pos, c, p_not_complete)

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
        if self.surface == "toric":
            self.measure_stabilizer_bulk(pos, c)
        elif self.surface == "planar":
            # Separate all stabilizers in bulk and boundaries.
            bulk_stabs = pos[:, self.plane[pos[0], pos[1]] == "o"]
            self.measure_stabilizer_bulk(bulk_stabs, c)
            self.measure_stabilizer_boundary(pos, c)

    def measure_stabilizer_bulk(self, pos, c):
            stab_qubits = self._stabilizer_qubits_bulk(pos)
            # Get all values on a multi dimensional array
            vals = self.qubits[c, stab_qubits[:, 0], stab_qubits[:, 1]]
            # Product over the desired dimension
            vals = np.prod(vals, axis=0)
            # Set the measurement results to the stabilizers
            self.qubits[0, pos[0], pos[1]] = vals

    def measure_stabilizer_boundary(self, pos, c):
            borders = ["t", "b", "l", "r"]
            for b in borders:
                self.measure_stabilizer_side(pos, b, c)

    def measure_stabilizer_side(self, pos, bord, c):
            # Separate all stabilizers in top, bottom, etc.
            bord_stabs = pos[:, self.plane[pos[0], pos[1]] == bord]

            # Get corresponding qubits
            bord_qubits = self._stabilizer_qubits_boundary(bord_stabs, bord)

            # Get all values on a multi dimensional array
            vals = self.qubits[c, bord_qubits[:, 0], bord_qubits[:, 1]]
            # Product over the desired dimension
            vals = np.prod(vals, axis=0)
            # Set the measurement results to the stabilizers
            self.qubits[0, bord_stabs[0], bord_stabs[1]] = vals


    def _stabilizer_qubits_bulk(self, pos):
        """Find qubits corresponing to the given stabilizers."""
        top = pos + np.array([[-1], [0]])
        bottom = pos + np.array([[1], [0]])
        left = pos + np.array([[0], [-1]])
        right = pos + np.array([[0], [1]])

        if self.surface == "toric":
            # Take the mod to account for cyclic boundaries
            # Top and left are automatically accounted for
            # when -1 is the index
            bottom = bottom % self.side
            right = right % self.side

        stab_qubits = np.stack((top, bottom, left, right), 0)
        return stab_qubits

    def _stabilizer_qubits_boundary(self, pos, bound):
        """Find qubits corresponing to the given stabilizers."""
        if bound == "t":
            a = pos + np.array([[1], [0]])
            b = pos + np.array([[0], [-1]])
            c = pos + np.array([[0], [1]])
        elif bound == "b":
            a = pos + np.array([[-1], [0]])
            b = pos + np.array([[0], [-1]])
            c = pos + np.array([[0], [1]])
        elif bound == "l":
            a = pos + np.array([[-1], [0]])
            b = pos + np.array([[1], [0]])
            c = pos + np.array([[0], [1]])
        elif bound == "r":
            a = pos + np.array([[-1], [0]])
            b = pos + np.array([[1], [0]])
            c = pos + np.array([[0], [-1]])

        stab_qubits = np.stack((a, b, c), 0)
        return stab_qubits

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
        self.qubits[0, self.tags == tag] *= lie

    def apply_measurement_error(self, p_lie):
        """Add measurement error to BOTH stars and plaqs stabilizers."""
        # Add measurement error to BOTH stabilizers
        self._stabilizer_lie("S", p_lie)
        self._stabilizer_lie("P", p_lie)

    def apply_qubit_error(self, pX, pZ):
        """Apply random error to the data qubits."""
        # Create the noise elements
        p = np.array([[pX], [pZ]])
        noise = 2*(np.random.rand(2, self.number_data_qubits) > p) - 1
        # Apply the noise
        self.qubits[:, self.tags == "Q"] *= noise

    def environmental_noise(self, p):
        pa = np.array([[p], [p]])/3.
        noise = 2*(np.random.rand(2, self.number_data_qubits) > pa) - 1
        # Apply the X Z noise
        self.qubits[:, self.tags == "Q"] *= noise

        # Apply the Y noise
        noise = 2*(np.random.rand(self.number_data_qubits) > p/3.) - 1
        self.qubits[0, self.tags == "Q"] *= noise
        self.qubits[1, self.tags == "Q"] *= noise


    def noisy_measurement_cycle(self, p_env):
        # Star measurements
        self.environmental_noise(p_env)
        self.noisy_measurement_specific(self.stars_round1, 0, "star")
        self.environmental_noise(p_env)
        self.noisy_measurement_specific(self.stars_round2, 0, "star")

        # Plaq measurements
        self.environmental_noise(p_env)
        self.noisy_measurement_specific(self.plaqs_round1, 1, "plaq")
        self.environmental_noise(p_env)
        self.noisy_measurement_specific(self.plaqs_round2, 1, "plaq")


    # NOTE: OLD function
    # def _apply_operation_error(self, pos, error):
    #     """Apply operation error."""
    #     pos_qubit1, pos_qubit2 = self.two_rand_stab_qubits(pos)
    #
    #     err_measurement, err_qubit1, err_qubit2 = error
    #     # Errors have shape:
    #     # [[e1X, e2X, e3X, ...], [e1Z, e2Z, e3Z, ...]]
    #     self.qubits[:, pos_qubit1[0], pos_qubit1[1]] *= err_qubit1
    #     self.qubits[:, pos_qubit2[0], pos_qubit2[1]] *= err_qubit2
    #
    #     # Apply error to stabilizer
    #     self.qubitts[0, pos[0], pos[1]] *= err_measurement

    def separate_bulk_boundary(self, pos):
        # Mask to identify bulk from boundary
        maskx = np.ones_like(pos[0], dtype=bool)
        masky = np.ones_like(pos[0], dtype=bool)
        # Top
        maskx[pos[0] == 0] = False
        # Bottom
        maskx[pos[0] == self.side - 1] = False
        # Left
        masky[pos[1] == 0] = False
        # Right
        masky[pos[1] == self.side - 1] = False

        bulk = np.concatenate((pos[:, maskx], pos[:, masky]), 1)
        boundary = np.concatenate((pos[:, np.invert(maskx)],
                                   pos[:, np.invert(masky)]), 1)
        return bulk, boundary

    def noisy_measurement(self, stabilizer):
        if stabilizer == "star":
            pos = self.stars
            c = 0
        elif stabilizer == "plaq":
            pos = self.plaqs
            c = 1

        self.noisy_measurement_specific(pos, c, stabilizer)

    def noisy_measurement_specific(self, pos, c, stabilizer):
        if self.surface == "toric":
            N = len(pos[0])
            m_err, q_err = self.errors.get_errors(N, stabilizer)
            stab_qubits = self._stabilizer_qubits_bulk(pos)
            # Apply error to qubits
            self.qubits[:, stab_qubits[:, 0], stab_qubits[:, 1]] *= q_err
            # Measure stabilizers
            self.measure_stabilizer(pos, c)
            # Apply errors to measurements
            self.qubits[0, pos[0], pos[1]] *= m_err
        elif self.surface == "planar":
            # First the bulk stabilizers
            bulk_stabs = pos[:, self.plane[pos[0], pos[1]] == "o"]
            N_bulk = len(bulk_stabs[0])
            m_err, q_err = self.errors.get_errors(N_bulk, stabilizer)
            bulk_qubits = self._stabilizer_qubits_bulk(bulk_stabs)

            # Apply error to qubits
            self.qubits[:, bulk_qubits[:, 0], bulk_qubits[:, 1]] *= q_err
            # Measure stabilizers
            self.measure_stabilizer(bulk_stabs, c)
            # Error to measurements
            self.qubits[0, bulk_stabs[0], bulk_stabs[1]] *= m_err

            # Now the boundaries
            for b in ["t", "b", "l", "r"]:
                bord_stabs = pos[:, self.plane[pos[0], pos[1]] == b]
                N_bord = len(bord_stabs[0])
                m_err, q_err = self.errors.get_errors(N_bord, stabilizer,
                                                      border=True)
                bord_qubits = self._stabilizer_qubits_boundary(bord_stabs, b)

                self.qubits[:, bord_qubits[:, 0], bord_qubits[:, 1]] *= q_err

                # Do measurement over border quibits - copied from function
                vals = self.qubits[c, bord_qubits[:, 0], bord_qubits[:, 1]]
                vals = np.prod(vals, axis=0)
                self.qubits[0, bord_stabs[0], bord_stabs[1]] = vals

                self.qubits[0, bord_stabs[0], bord_stabs[1]] *= m_err

    def _select_stabilizer(self, stabilizer):
        """Return useful parameters for stabilizer type."""
        if stabilizer == "star":
            c = 0
        if stabilizer == "plaq":
            c = 1

        return c

    # def plot(self, stabilizer, backup=False):
    #     """Plot the surface code."""
    #     if stabilizer == "star":
    #         data = self.qubits[0].copy()
    #         data[self.tags == "S"] *= 2
    #         data[self.tags == "P"] = 1
    #     if stabilizer == "plaq":
    #         data = self.qubits[0].copy()
    #         data[self.tags == "Q"] = self.qubits[1, self.tags == "Q"]
    #         data[self.tags == "P"] *= 2
    #         data[self.tags == "S"] = 1
    #
    #     # Return data to plot
    #     # return data, self.cmap, self.cmap_norm
    #     plt.figure()
    #     plt.imshow(data, cmap=self.cmap, norm=self.cmap_norm)
    #     # plt.colorbar()
    #     # plt.show()

    def reset(self):
        """Reset surface code to default configuration."""
        self.qubits.fill(1)
        self.qubits[1, self.tags != "Q"] = 0

    def get_stars(self):
        """Get the values of  all star qubits."""
        return self.qubits[0, self.tags == "S"]

    def get_plaqs(self):
        """Get the values of all plaq qubits."""
        return self.qubits[0, self.tags == "P"]

    def measure_logical(self):
        """Meausure the logical qubits."""
        if self.surface == "toric":
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


    def correct_error(self, error_type, match, time=0):
        if len(match) == 0:
            return

        c = self._select_stabilizer(error_type)
        m = self.side

        # Index where one pair is on the last time sheet and the other
        # in virtual time
        # ind = np.invert(np.prod(match[:, :, 2] == time, 1)).astype(bool)
        # faulty_stabs = match[ind, 0, :].transpose()
        # print("Faulty Stabs")
        # print(faulty_stabs)
        # self.qubits[0, faulty_stabs[0], faulty_stabs[1]] *= -1
        #
        # match = match[np.invert(ind)]

        if self.surface == "planar":
            for pair in match:
                px, py, pt = pair[0]
                qx, qy, qt = pair[1]

                # Fix faulty measurements
                # if pt == time + 1 or qt == time + 1:
                #     if abs(pt - qt) == 2 and px == qx and py == qy:
                #         # print("Faulty measurement: ", px, py)
                #         self.qubits[0, px, py] *= -1
                #         continue

                dx = qx - px
                sx = np.sign(dx)
                dx = np.abs(dx)
                dy = qy - py
                sy = np.sign(dy)
                dy = np.abs(dy)

                stepsx = sx*np.arange(1, dx, 2) + px
                stepsy = sy*np.arange(1, dy, 2) + py
                coord_x = np.ones_like(stepsy) * qx
                coord_y = np.ones_like(stepsx) * py
                stepsx = np.append(stepsx, coord_x)
                stepsy = np.append(coord_y, stepsy)

                # Apply error correction path
                # print("STEPS")
                # print(stepsx)
                # print(stepsy)
                self.qubits[c, stepsx, stepsy] *= -1

        elif self.surface == "toric":
            for pair in match:
                # print("Pair:", pair)
                px, py, pt = pair[0]
                qx, qy, qt = pair[1]

                # Fix faulty measurements
                # if pt == time + 1 or qt == time + 1:
                #     if abs(pt - qt) == 2 and px == qx and py == qy:
                #         # print("Faulty measurement: ", px, py)
                #         self.qubits[0, px, py] *= -1
                #         continue

                dx = (qx - px) % m
                dy = (qy - py) % m

                if dx < self.distance:
                    endx = qx
                    stepsx = np.arange(1, dx, 2)
                    stepsx = (stepsx + px) % m
                    coord_y = np.ones_like(stepsx) * py
                else:
                    endx = px
                    stepsx = np.arange(1, m - dx, 2)
                    stepsx = (stepsx + qx) % m
                    coord_y = np.ones_like(stepsx) * qy

                if dy < self.distance:
                    stepsy = (np.arange(1, dy, 2) + py) % m

                else:
                    stepsy = (np.arange(1, m - dy, 2) + qy) % m

                coord_x = np.ones_like(stepsy) * endx

                stepsx = np.append(stepsx, coord_x).astype(int)
                stepsy = np.append(coord_y, stepsy).astype(int)

                # print("Steps")
                # print(stepsx)
                # print(stepsy)
                self.qubits[c, stepsx, stepsy] *= -1
