"""
Complete list of errors for the noisy operators.
"""
I, X, Y, Z = [1, 1], [-1, 1], [-1, -1], [1, -1]

errorPlaquete1 = [[1, [I, I]], [1, [I, Z]], [1, [I, X]], [1, [I, Y]],
                   [1, [X, X]], [1, [X, Y]], [1, [Y, Y]], [1, [X, Z]],
                  [1, [Y, Z]], [1, [Z, Z]], [-1,[I,I]],[-1,[I,Z]],
                  [-1,[I,X]],[-1,[I,Y]],[-1,[X,X]],[-1,[X,Y]],
                  [-1,[Y,Y]],[-1,[X,Z]],[-1,[Y,Z]],[-1,[Z,Z]]]

errorPlaquete2 = [[1,[I,I]],[1,[I,Z]],[-1,[I,X]],[-1,[I,Y]],
                  [1,[X,X]],[1,[X,Y]],[1,[Y,Y]],[-1,[X,Z]],
                  [-1,[Y,Z]],[1,[Z,Z]],[-1,[I,I]],[-1,[I,Z]],
                  [1,[I,X]],[1,[I,Y]],[-1,[X,X]],[-1,[X,Y]],
                  [-1,[Y,Y]],[1,[X,Z]],[1,[Y,Z]],[-1,[Z,Z]]]

errorStar1 = np.array([[1,[I,I]],[1,[I,X]],[1,[I,Z]],[1,[I,Y]],
              [1,[Z,Z]],[1,[Z,Y]],[1,[Y,Y]],[1,[Z,X]],
              [1,[Y,X]],[1,[X,X]],[-1,[I,I]],[-1,[I,X]],
              [-1,[I,Z]],[-1,[I,Y]],[-1,[Z,Z]],[-1,[Z,Y]],
              [-1,[Y,Y]],[-1,[Z,X]],[-1,[Y,X]],[-1,[X,X]]])

errorStar2 = [[1,[I,I]],[1,[I,X]],[-1,[I,Z]],[-1,[I,Y]],
            [1,[Z,Z]],[1,[Z,Y]],[1,[Y,Y]],[-1,[Z,X]],
            [-1,[Y,X]],[1,[X,X]],[-1,[I,I]],[-1,[I,X]],
            [1,[I,Z]],[1,[I,Y]],[-1,[Z,Z]],[-1,[Z,Y]],
            [-1,[Y,Y]],[1,[Z,X]],[1,[Y,X]],[-1,[X,X]]]



def processErrors(errorProbabilities, errorList):
            sortList = sorted(zip(errorProbabilities, errorList),reverse=True)
            # Separate errors and probabilities
            probs, errors = zip(*sortList)
            # Cumulative probalities of the errors
            cumulativeProbs = np.cumsum(np.array(probs1))
