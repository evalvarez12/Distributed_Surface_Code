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


errTestVec = [0.926639, 0.006904, 0.003904, 0.003904, 0.000012, 0.000024, 0.000012,
              0.000048, 0.000048, 0.001062, 0.042477, 0.006868, 0.003904, 0.003904,
              0.000012, 0.000024, 0.000012, 0.000048, 0.000048, 0.000147]

def processErrors(errorProbabilities, errorList):
            # TODO whats the point of this?
            sortList = sorted(zip(errorProbabilities, errorList), reverse=True)
            # Separate errors and probabilities
            probs, errors = zip(*sortList)
            # Cumulative probalities of the errors
            cumulativeProbs = np.cumsum(np.array(probs))

            numberErrors = len(errorList)
            measurementErr = np.zeros(numberErrors)
            errQubit1 = np.zeros((numberErrors,2))
            errQubit2 = np.zeros((numberErrors,2))
            for i in range(numberErrors):
                measurementErr[i] = errorList[i][0]
                errQubit1[i] = errorList[i][1][0]
                errQubit2[i] = errorList[i][1][1]

            # Errors in shape:
            # [[e1X, e2X, e3X, ...], [e1Z, e2Z, e3Z, ...]]
            errQubit1 = errQubit1.transpose()
            errQubit2 = errQubit2.transpose()
            errors = (measurementErr, errQubit1, errQubit2)
            return cumulativeProbs, errors
