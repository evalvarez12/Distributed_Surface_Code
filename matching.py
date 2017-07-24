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

    if anyons  == 0:
        return []

    graph = makeGraph(anyons, weights)

    numberNodes = len(anyons)
    print(graph)
    matching = pm.getMatching(numberNodes, graph)
    return(matching)
    # matching_pairs=[[i,matching[i]] for i in range(numberNodesnodes) if matching[i]>i]
    # points=[] if len(matching_pairs)==0 else [[nodes_list[i] for i in x] for x in matching_pairs]
    #
    # return points


def makeGraph(nodes, weights=[1, 1]):
    # TODO redo this
    numberNodes = len(nodes)
    m = 2*numberNodes
    wT, wS = weights
    graph = []
    for i in range(numberNodes-1):
        p1, p2, p3 = nodes[i]
        for j in range(i, numberNodes-1):
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
