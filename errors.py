"""
Complete list of processed errors for the noisy operators.
"""
import numpy as np
import pickle
from os.path import dirname, realpath
# import decomposition.generate as gen

I, X, Y, Z = [1, 1], [-1, 1], [-1, -1], [1, -1]

class Generator:
    def __init__(self, surface, ps, pm, pg, eta, a0, a1, theta, protocol):
        # self.generator = gen.Generator()

        chi_star = self._load_model(ps, pm, pg, eta, a0, a1, theta,
                                    4, "X", protocol)
        chi_plaq = self._load_model(ps, pm, pg, eta, a0, a1, theta,
                                    4, "Z", protocol)
        self.chi = [chi_star, chi_plaq]

        self.chi_keys = [np.array(list(chi_star.keys())),
                         np.array(list(chi_plaq.keys()))]
        self.chi_vals = [np.array(list(chi_star.values())),
                         np.array(list(chi_plaq.values()))]
        self.errors = [self._symbol_to_error_list(self.chi_keys[0]),
                       self._symbol_to_error_list(self.chi_keys[1])]
        self.indexes = range(len(self.chi_keys[0]))

        if surface == "planar":
            chi_star_border = self._load_model(ps, pm, pg, eta, a0, a1, theta,
                                               3, "X", protocol)
            chi_plaq_border = self._load_model(ps, pm, pg, eta, a0, a1, theta,
                                               3, "Z", protocol)
            self.chi_border = [chi_star_border, chi_plaq_border]

            self.chi_keys_border = [np.array(list(chi_star_border.keys())),
                                    np.array(list(chi_plaq_border.keys()))]
            self.chi_vals_border = [np.array(list(chi_star_border.values())),
                                    np.array(list(chi_plaq_border.values()))]
            self.errors_border = [self._symbol_to_error_list(self.chi_keys_border[0], True),
                                  self._symbol_to_error_list(self.chi_keys_border[1], True)]
            self.indexes_border = range(len(self.chi_keys_border[0]))

    def _generate_name(self, ps, pm, pg, eta, a0, a1, theta, stab_size, parity, protocol):
        """Name for CHI matrix file."""
        param_names = ["ps=" + str(ps), "pm=" + str(pm),
                       "pg=" + str(pg), "eta=" + str(round(eta, 4)),
                       "a0=" + str(round(a0, 4)), "a1=" + str(round(a1, 4)),
                       "theta=" + str(theta)]

        param_names = "_".join(param_names)
        file_name = ["CHI", protocol, parity, str(stab_size)]
        file_name = "_".join(file_name)
        # The address of the parent directory
        script_path = dirname(dirname(realpath(__file__)))
        file_name = (script_path + "/data/" + file_name
                     + "_" + param_names + ".dict")
        return file_name

    def _load_model(self, ps, pm, pg, eta, a0, a1, theta, stab_size, parity, protocol):
        file_name = self._generate_name(ps, pm, pg, eta, a0, a1, theta,
                                        stab_size, parity, protocol)
        # try:
        pickle_in = open(file_name, "rb")
        return pickle.load(pickle_in)
        # except:
            # raise NameError("Model file not found: " + file_name)

    def get_errors(self, num_errors, stabilizer, border=False):
        if stabilizer == "star":
            c = 0
        elif stabilizer == "plaq":
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

    def _symbol_to_error_list(self, symbol_list, border=False):
        N = len(symbol_list)
        if border:
            d_qubits = 3
        else:
            d_qubits = 4

        qubit_errs = np.ones((d_qubits, 2, N))
        measurement_errs = np.ones(N)
        for i in range(N):
            m, q = self._symbol_to_error(symbol_list[i])

            # NOTE: Shuffle to take a random permutation of the error
            np.random.shuffle(q)

            measurement_errs[i] *= m
            qubit_errs[:, :, i] *= q

        return [measurement_errs, qubit_errs]

    def _symbol_to_error(self, symbol):
        if "N" in symbol:
            measurement = -1
            symbol = symbol.replace("_NOK", "")
        else:
            measurement = 1
            symbol = symbol.replace("_OK", "")

        l_symbol = list(symbol)
        errors = [self._pauli_error(s) for s in l_symbol]
        return measurement, errors

    def _pauli_error(self, s):
        if s == "I":
            return I
        elif s == "X":
            return X
        elif s == "Y":
            return Y
        elif s == "Z":
            return Z
