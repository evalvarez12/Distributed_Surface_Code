"""
Complete list of errors for the noisy operators.
"""
import numpy as np
import decomposition.generate as gen

I, X, Y, Z = [1, 1], [-1, 1], [-1, -1], [1, -1]

class Generator:
    def __init__(self, surface, ps, pm, pg, pn, protocol):
        self.generator = gen.Generator()

        chi_star = self.generator.ask_model(ps, pm, pg, pn, 4,
                                            "X", protocol)
        chi_plaq = self.generator.ask_model(ps, pm, pg, pn, 4,
                                            "Z", protocol)
        # self.chi = [chi_star, chi_plaq]

        self.chi_keys = [np.array(list(chi_star.keys())),
                         np.array(list(chi_plaq.keys()))]
        self.chi_vals = [np.array(list(chi_star.values())),
                         np.array(list(chi_plaq.values()))]
        self.errors = [self.symbol_to_error_list(self.chi_keys[0]),
                       self.symbol_to_error_list(self.chi_keys[1])]
        self.indexes = range(len(self.chi_keys[0]))

        if surface == "planar":
            chi_star_border = self.generator.ask_model(ps, pm, pg, pn, 3,
                                                       "X", protocol)
            chi_plaq_border = self.generator.ask_model(ps, pm, pg, pn, 3,
                                                       "Z", protocol)
            # self.chi_border = [chi_star_border, chi_plaq_border]

            self.chi_keys_border = [np.array(list(chi_star_border.keys())),
                                    np.array(list(chi_plaq_border.keys()))]
            self.chi_vals_border = [np.array(list(chi_star_border.values())),
                                    np.array(list(chi_plaq_border.values()))]
            self.errors_border = [self.symbol_to_error_list(self.chi_keys_border[0], True),
                                  self.symbol_to_error_list(self.chi_keys_border[1], True)]
            self.indexes_border = range(len(self.chi_keys_border[0]))

    def get_errors(self, num_errors, stabilizer, border=False):
        if stabilizer == "X":
            c = 0
        elif stabilizer == "Z":
            c = 1

        if border:
            i = self.indexes_border
            v = self.chi_vals_border[c]
            e = self.errors_border[c]
        else:
            i = self.indexes
            v = self.chi_vals[c]
            e = self.errors[c]

        e_index = np.array(np.random.choice(i, num_errors, p=v))
        m_errors = e[0][e_index]
        # Qubit erros format:
        # array([[[ 1., -1.],
        #         [ 1., -1.]],
        #
        #        [[ 1.,  1.],
        #         [ 1.,  1.]],   All 4 data qubits for each stabilizer
        #                        q_errors[0] gives the error for all qubits
        #        [[ 1., -1.],    on top of the stabilizer.
        #         [ 1.,  1.]],   q_errors[1] for the bottom and so on.
        #
        #        [[ 1.,  1.],
        #         [ 1.,  1.]]])
        q_errors = e[1][:, :, e_index]

        # Swap axes to make it manageable
        q_errors = np.swapaxes(q_errors, 0, 1)
        return m_errors, q_errors

    def symbol_to_error_list(self, symbol_list, border=False):
        N = len(symbol_list)
        if border:
            d_qubits = 3
        else:
            d_qubits = 4

        qubit_errs = np.ones((d_qubits, 2, N))
        measurement_errs = np.ones(N)
        for i in range(N):
            m, q = self.symbol_to_error(symbol_list[i])

            # NOTE: Shuffle to take a random permutation of the error
            np.random.shuffle(q)

            measurement_errs[i] *= m
            qubit_errs[:, :, i] *= q

        return [measurement_errs, qubit_errs]

    def symbol_to_error(self, symbol):
        if "N" in symbol:
            measurement = -1
            symbol = symbol.replace("_NOK", "")
        else:
            measurement = 1
            symbol = symbol.replace("_OK", "")

        l_symbol = list(symbol)
        errors = [self.pauli_error(s) for s in l_symbol]
        return measurement, errors

    def pauli_error(self, s):
        if s == "I":
            return I
        elif s == "X":
            return X
        elif s == "Y":
            return Y
        elif s == "Z":
            return Z






#
#
# errorPlaquete1 = [[1, [I, I]], [1, [I, Z]], [1, [I, X]], [1, [I, Y]],
#                   [1, [X, X]], [1, [X, Y]], [1, [Y, Y]], [1, [X, Z]],
#                   [1, [Y, Z]], [1, [Z, Z]], [-1, [I, I]], [-1, [I, Z]],
#                   [-1, [I, X]], [-1, [I, Y]], [-1, [X, X]], [-1, [X, Y]],
#                   [-1, [Y, Y]], [-1, [X, Z]], [-1, [Y, Z]], [-1, [Z, Z]]]
#
# errorPlaquete2 = [[1, [I, I]], [1, [I, Z]], [-1, [I, X]], [-1, [I, Y]],
#                   [1, [X, X]], [1, [X, Y]], [1, [Y, Y]], [-1, [X, Z]],
#                   [-1, [Y, Z]], [1, [Z, Z]], [-1, [I, I]], [-1, [I, Z]],
#                   [1, [I, X]], [1, [I, Y]], [-1, [X, X]], [-1, [X, Y]],
#                   [-1, [Y, Y]], [1, [X, Z]], [1, [Y, Z]], [-1, [Z, Z]]]
#
# errorStar1 = [[1, [I, I]], [1, [I, X]], [1, [I, Z]], [1, [I, Y]],
#               [1, [Z, Z]], [1, [Z, Y]], [1, [Y, Y]], [1, [Z, X]],
#               [1, [Y, X]], [1, [X, X]], [-1, [I, I]], [-1, [I, X]],
#               [-1, [I, Z]], [-1, [I, Y]], [-1, [Z, Z]], [-1, [Z, Y]],
#               [-1, [Y, Y]], [-1, [Z, X]], [-1, [Y, X]], [-1, [X, X]]]
#
# errorStar2 = [[1, [I, I]], [1, [I, X]], [-1, [I, Z]], [-1, [I, Y]],
#               [1, [Z, Z]], [1, [Z, Y]], [1, [Y, Y]], [-1, [Z, X]],
#               [-1, [Y, X]], [1, [X, X]], [-1, [I, I]], [-1, [I, X]],
#               [1, [I, Z]], [1, [I, Y]], [-1, [Z, Z]], [-1, [Z, Y]],
#               [-1, [Y, Y]], [1, [Z, X]], [1, [Y, X]], [-1, [X, X]]]
#
#
# err_test_vec = [0.926639, 0.006904, 0.003904, 0.003904, 0.000012, 0.000024, 0.000012,
#               0.000048, 0.000048, 0.001062, 0.042477, 0.006868, 0.003904, 0.003904,
#               0.000012, 0.000024, 0.000012, 0.000048, 0.000048, 0.000147]
#
# def process_errors(errorProbabilities, errorList):
#             # TODO whats the point of this?
#             sort_list = sorted(zip(error_probabilities, error_list), reverse=True)
#             # Separate errors and probabilities
#             probs, errors = zip(*sort_list)
#             # Cumulative probalities of the errors
#             cumulative_probs = np.cumsum(np.array(probs))
#
#             number_errors = len(error_list)
#             measurement_err = np.zeros(number_errors)
#             err_qubit1 = np.zeros((number_errors,2))
#             err_qubit2 = np.zeros((number_errors,2))
#             for i in range(number_errors):
#                 measurement_err[i] = error_list[i][0]
#                 err_qubit1[i] = error_list[i][1][0]
#                 err_qubit2[i] = error_list[i][1][1]
#
#             # Errors in shape:
#             # [[e1X, e2X, e3X, ...], [e1Z, e2Z, e3Z, ...]]
#             err_qubit1 = err_qubit1.transpose()
#             err_qubit2 = err_qubit2.transpose()
#             errors = (measurement_err, err_qubit1, err_qubit2)
#             return cumulative_probs, errors
