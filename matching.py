import blossom5.pyMatch as pm
import numpy as np

def match_toric3D(size, anyons, weights=[1,1]):
    # TODO consideration when toroid and when planar
    if len(anyons)  == 0:
        return []

    graph = make_graph_toric(size, anyons, weights)
    # print(graph)
    # matching: indexes to which anyon conects
    number_nodes = len(anyons)
    matching = pm.getMatching(number_nodes, graph)
    # print(matching)
    pairs_ind = [[i, matching[i]] for i in range(number_nodes) if matching[i]>i]
    # print(pairs_ind)
    pairs = [] if len(pairs_ind)==0 else [anyons[p] for p in pairs_ind]
    # Pairs format:
    # np.array([[pair1, pair2], [pair1, pair2], ...])
    return np.array(pairs)


def make_graph_toric(size, nodes, weights=[1, 1]):
    # Nodes array format is:
    # [[x1, y1, t1], [x2, y2, t2], [x3, y3, t3], ...]
    number_nodes = len(nodes)
    # m is used to find shortest distance accros the toroid
    m = 2*size + 1

    # Spatial and time weights
    wt, ws = weights

    graph = []
    # TODO do without this loops?
    for i in range(number_nodes-1):
        px, py, pt = nodes[i]
        for j in range(i+1, number_nodes):
            qx, qy, qt = nodes[j]

            difft = (qt - pt)*wt
            diffx = abs(qx - px)
            diffx = min([diffx, m-diffx])*ws
            diffy = abs(qy - py)
            diffy = min([diffy, m-diffy])*ws

            weight = difft + diffx + diffy
            graph += [[i, j, weight]]
    # List of values for time distance weights
    return graph
