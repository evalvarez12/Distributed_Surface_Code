"""
Load a error model for the specified parameters.
The error model must be in the folder data/ following the name format
stablished by decomposition/tools/names.py.

author: Eduardo Villasenor
created-on: 21/08/17
"""
import numpy as np
import pickle
from os.path import dirname, realpath
import decomposition.tools.names as names

# Pauli errors on stabilizer representation
I, X, Y, Z = [1, 1], [1, -1], [-1, -1], [-1, 1]


class Generator:
    def __init__(self, surface, ps, pm, pg, eta, a0, a1, theta, protocol):
        """Init function.

        Parameters
        ----------
        surface : (string) planar or toric code
        ps : (scalar) single qubit gate error rate.
        pm : (scalar) measurement error rate.
        pg : (scalar) two qubit gate error rate.
        eta : (scalar) detection efficiency.
        a0 : (scalar) extra environmental error when electron spin is being operated.
        a1 : (scalar) default environmental error.
        theta : (scalar) determines how the states are initialized when generating remote
                entanglement.
        protocol : (string) name of the protocol used in generating the states
        """

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
        file_name = names.chi(ps, pm, pg, eta, a0, a1, theta,
                              stab_size, parity, protocol)
        return file_name

    def _load_model(self, ps, pm, pg, eta, a0, a1, theta, stab_size, parity, protocol):
        # Load a file with a error model in a dictionary
        file_name = self._generate_name(ps, pm, pg, eta, a0, a1, theta,
                                        stab_size, parity, protocol)
        # try:
        pickle_in = open(file_name, "rb")
        return pickle.load(pickle_in)
        # except:
            # raise NameError("Model file not found: " + file_name)

    def get_errors(self, num_errors, stabilizer, border=False):
        """
        Get a list of errors corresponding to the loaded error model.

        Paramaters
        -----------
        num_errors : (int) number of stabilizers measured on which errors are
                     applied
        stabilizer : (string) plaq or star, type of the stabilizer measured
        border : (bool) if the stabilizer measured is in the border of the
                 planar code
        """
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

        # Shuffle qubit errors to restore the lost permutations
        q_errors = self._shuffle_qubit_error(q_errors)

        # Swap axes to make it manageable
        q_errors = np.swapaxes(q_errors, 0, 1)
        return m_errors, q_errors

    def _shuffle_qubit_error(self, err_list):
        # NOTE: shuffle acts on pointers so this function can be inserted in
        # code for more efficiency
        # Change shape to shuffle correct axis
        err_list = err_list.transpose((1, 0, 2))
        # Shuffle qubit errors to restore the lost permutations
        np.random.shuffle(err_list)
        # Recover shape
        err_list = err_list.transpose((1, 0, 2))
        return err_list

    def _symbol_to_error_list(self, symbol_list, border=False):
        # Transform the symbols in the error dictionary to lists with the
        # corresponding errors
        N = len(symbol_list)
        if border:
            d_qubits = 3
        else:
            d_qubits = 4

        # Create lists to save the errors
        qubit_errs = np.ones((d_qubits, 2, N))
        measurement_errs = np.ones(N)
        # Iterate through the error filling the lists with the corresponding
        # values
        for i in range(N):
            m, q = self._symbol_to_error(symbol_list[i])

            measurement_errs[i] *= m
            qubit_errs[:, :, i] *= q

        return [measurement_errs, qubit_errs]

    def _symbol_to_error(self, symbol):
        # Tansform a errror symbol into a list representation
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
        # Individual string to list representation
        if s == "I":
            return I
        elif s == "X":
            return X
        elif s == "Y":
            return Y
        elif s == "Z":
            return Z
