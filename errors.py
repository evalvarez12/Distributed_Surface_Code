"""
Complete list of errors for the noisy operators.
"""
import numpy as np
import decomposition.generate

I, X, Y, Z = [1, 1], [-1, 1], [-1, -1], [1, -1]

class Errors:
    def __init__(self, surface, ps, pm, pg, pn, protocol):
        self.generator = generate.Generator()

        chi_star = self.generator.ask_model(ps, pm, pg, pn, 4,
                                            "X", protocol)
        chi_plaq = self.generator.ask_model(ps, pm, pg, pn, 4,
                                            "Z", protocol)
        self.chi = [chi_star, chi_plaq]
        self.errors = {}

        # for k in self.chi.values():
        #     self.errors[k] = self.symbol_to_error(k)


        if surface == "planar":
            chi_star_border = self.generator.ask_model(ps, pm, pg, pn, 4,
                                                       "X", protocol)
            chi_plaq_border = self.generator.ask_model(ps, pm, pg, pn, 4,
                                                       "Z", protocol)
            self.chi_border = [chi_star_border, chi_plaq_border]
            self.error_border = {}

            # for k in self.chi_border.values():
            #     self.errors_border[k] = self.symbol_to_error(k)

def process_chi(self):
    for k in self.chi.keys():
        self.error[k] = symbol_to_error(k)

def symbol_to_error(self, symbol):
    measurement = 1
    if "N" in symbol:
        measurement = -1

    symbol = symbol.replace("_NOK", "")
    l_symbol = list(symbol)
    errors = [pauli_error(s) for s in l_symbol]
    return measurement, errors

def get_errors(self, stab_positions, stabilizer, border=False):
    errors_sym = self.choose_error(len(stab_positions), stabilizer, border)
    if border:
        no_err = "III_OK"
        d_qubits = 3
    else:
        no_err = "IIII_OK"
        d_qubits = 4

    # Remove where no errors ocurred
    stab_positions = stab_positions[errors_sym != no_err]
    errors_sym = errors_sym[errors_sym != no_err]

    # Arrays to save errors
    N = len(errors_sym)
    measurement_errors = np.ones(N)
    qubit_errors = np.ones((d_qubits, 2, N))
    qubit1
    for i in range(N):
        m_err, q_err = self.symbol_to_error(errors_sym[i])
        measurement_errors[i] *= m_err
        # NOTE: Shuffle the errors on the data qubits
        # np.random.shuffle(q_err)

        # Get all the errors of the data qubits for a given stabilizer
        qubit_errors[:, :, i] *= q_err

    return stab_positions, measurement_errors, qubit_errors

def choose_error(self, num_errors, stabilizer, border):
    if stabilizer == "X":
        c = 0
    elif stabilizer == "Z":
        c = 1

    if border:
        k = list(self.chi_border[c].keys())
        v = list(self.chi_border[c].values())
    else:
        k = list(self.chi[c].keys())
        v = list(self.chi[c].values())

    errors_keys = np.random.choice(k, num_errors, p=v)
    return errors_keys

def pauli_error(s):
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
