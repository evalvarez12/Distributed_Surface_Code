import blossom5.pyMatch as pm


# TODO is this important?
# global time_lookup, weight_lookup
#
# timeCutOff = 15
# spaceCutOff = 2 * 20
#
# time_lookup = {}
# for t in range(cutoff + 1):
#     time_lookup[t] = t
#
# weight_lookup = {}
# for i in range(m):
#     weight_lookup[i] = {}
#     for j in range(m):
#         diff = abs(i - j)
#         weight_lookup[i][j] = min([diff, m-diff])



def matchToric3D(size, anyons, weights=[1,1]):

    if len(anyons)  == 0:
        return []

    graph = makeGraph(anyons, weights)
    print(anyons)
    numberNodes = len(anyons)
    print(graph)
    # matching: indexes to which anyon conects
    matching = pm.getMatching(numberNodes, graph)
    print(matching)
    pairsInd = [[i, matching[i]] for i in range(numberNodes) if matching[i]>i]
    print(pairsInd)
    pairs = [] if len(pairsInd)==0 else [[anyons[p]] for p in pairsInd]

    return pairs


def makeGraph(nodes, weights=[1, 1]):
    # TODO redo this
    numberNodes = len(nodes)
    m = 2*numberNodes
    wT, wS = weights
    graph = []

    # ind = np.arange(numberNodes)
    # graph = [nodes - nodes[i+1:] for i in range(numberNodes)]
    # graphIdexes = []
    for i in range(numberNodes-1):
        p1, p2, p3 = nodes[i]
        for j in range(i+1, numberNodes):
            q1, q2, q3 = nodes[j]

            wt = (q3 - p3)*wT
            diffx = abs(q1 - p1)
            diffx = min([diffx, m-diffx])*wS
            diffy = abs(q2 - q2)
            diffy = min([diffy, m-diffy])*wS

            weight = wt + diffx + diffy
            graph += [[i, j, weight]]
    # List of values for time distance weights
    return graph

def weightNorm(pos, m, weights):
    # TODO: use this?
    # pos = [[px, py, pt], [px, py, pt]]
    swT, wS = weights
    res = abs(pos)
    res = (pos[:, 0]%m + pos[:, 1]%m)*wS + pos[:, 2]*wT
    return res
