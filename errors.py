"""
Complete list of errors for the noisy operators.
"""
import numpy as np

I, X, Y, Z = [1, 1], [-1, 1], [-1, -1], [1, -1]

errorPlaquete1 = [[1, [I, I]], [1, [I, Z]], [1, [I, X]], [1, [I, Y]],
                  [1, [X, X]], [1, [X, Y]], [1, [Y, Y]], [1, [X, Z]],
                  [1, [Y, Z]], [1, [Z, Z]], [-1, [I, I]], [-1, [I, Z]],
                  [-1, [I, X]], [-1, [I, Y]], [-1, [X, X]], [-1, [X, Y]],
                  [-1, [Y, Y]], [-1, [X, Z]], [-1, [Y, Z]], [-1, [Z, Z]]]

errorPlaquete2 = [[1, [I, I]], [1, [I, Z]], [-1, [I, X]], [-1, [I, Y]],
                  [1, [X, X]], [1, [X, Y]], [1, [Y, Y]], [-1, [X, Z]],
                  [-1, [Y, Z]], [1, [Z, Z]], [-1, [I, I]], [-1, [I, Z]],
                  [1, [I, X]], [1, [I, Y]], [-1, [X, X]], [-1, [X, Y]],
                  [-1, [Y, Y]], [1, [X, Z]], [1, [Y, Z]], [-1, [Z, Z]]]

errorStar1 = [[1, [I, I]], [1, [I, X]], [1, [I, Z]], [1, [I, Y]],
              [1, [Z, Z]], [1, [Z, Y]], [1, [Y, Y]], [1, [Z, X]],
              [1, [Y, X]], [1, [X, X]], [-1, [I, I]], [-1, [I, X]],
              [-1, [I, Z]], [-1, [I, Y]], [-1, [Z, Z]], [-1, [Z, Y]],
              [-1, [Y, Y]], [-1, [Z, X]], [-1, [Y, X]], [-1, [X, X]]]

errorStar2 = [[1, [I, I]], [1, [I, X]], [-1, [I, Z]], [-1, [I, Y]],
              [1, [Z, Z]], [1, [Z, Y]], [1, [Y, Y]], [-1, [Z, X]],
              [-1, [Y, X]], [1, [X, X]], [-1, [I, I]], [-1, [I, X]],
              [1, [I, Z]], [1, [I, Y]], [-1, [Z, Z]], [-1, [Z, Y]],
              [-1, [Y, Y]], [1, [Z, X]], [1, [Y, X]], [-1, [X, X]]]


err_test_vec = [0.926639, 0.006904, 0.003904, 0.003904, 0.000012, 0.000024, 0.000012,
              0.000048, 0.000048, 0.001062, 0.042477, 0.006868, 0.003904, 0.003904,
              0.000012, 0.000024, 0.000012, 0.000048, 0.000048, 0.000147]

def process_errors(errorProbabilities, errorList):
            # TODO whats the point of this?
            sort_list = sorted(zip(error_probabilities, error_list), reverse=True)
            # Separate errors and probabilities
            probs, errors = zip(*sort_list)
            # Cumulative probalities of the errors
            cumulative_probs = np.cumsum(np.array(probs))

            number_errors = len(error_list)
            measurement_err = np.zeros(number_errors)
            err_qubit1 = np.zeros((number_errors,2))
            err_qubit2 = np.zeros((number_errors,2))
            for i in range(number_errors):
                measurement_err[i] = error_list[i][0]
                err_qubit1[i] = error_list[i][1][0]
                err_qubit2[i] = error_list[i][1][1]

            # Errors in shape:
            # [[e1X, e2X, e3X, ...], [e1Z, e2Z, e3Z, ...]]
            err_qubit1 = err_qubit1.transpose()
            err_qubit2 = err_qubit2.transpose()
            errors = (measurement_err, err_qubit1, err_qubit2)
            return cumulative_probs, errors
