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

    graph = makeGraphToric(size, anyons, weights)
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


def makeGraphToric(size, nodes, weights=[1, 1]):
    numberNodes = len(nodes)
    # m is used to find shortest distance accros the toroid
    m = 2*size + 1

    # Spatial and time weights
    wT, wS = weights

    graph = []

    # ind = np.arange(numberNodes)
    # graph = [nodes - nodes[i+1:] for i in range(numberNodes)]
    # graphIdexes = []
    for i in range(numberNodes-1):
        px, py, pt = nodes[i]
        for j in range(i+1, numberNodes):
            qx, qy, qt = nodes[j]

            difft = (qt - pt)*wT
            diffx = abs(qx - px)
            diffx = min([diffx, m-diffx])*wS
            diffy = abs(qy - qy)
            diffy = min([diffy, m-diffy])*wS

            weight = difft + diffx + diffy
            graph += [[i, j, weight]]
    # List of values for time distance weights
    return graph
