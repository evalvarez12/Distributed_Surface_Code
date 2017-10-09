import blossom5.pyMatch as pm
import numpy as np


def match_toric_3D(size, anyons, weights=[1, 1]):
    """
    Find a matching to fix the errors in a 3D planar code given the positions
    of '-1' stabilizer outcomes

    Parameters:
    -----------
    size -- The dimension of the code
    anyons -- A list of the locations of all '-1' value stabilizers.
              [[x0,y0,t0],[x1,y1,t1],...]
    weights -- The multiplicative weighting that should be assigned to graph
               edges in the [space,time] dimensions. Default: [1,1]

    Returns:
    --------
    A list containing all the input anyon positions grouped into pairs. [[[x0,y0,t0],[x1,y1,t1]],[[x2,y2,t2],...
    """
    # TODO consideration when toroid and when planar
    if len(anyons) == 0:
        return []

    graph = make_graph_toric(size, anyons, weights)
    # print(graph)
    # matching: indexes to which anyon conects
    number_nodes = len(anyons)
    if number_nodes % 2 == 1:
        raise ValueError("Number of nodes is odd!")
    matching = pm.getMatching(number_nodes, graph)
    # print(matching)
    pairs_ind = [[i, matching[i]] for i in range(number_nodes) if matching[i] > i]
    # print(pairs_ind)
    pairs = [] if len(pairs_ind) == 0 else [anyons[p] for p in pairs_ind]
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

def match_planar_3D(size, anyons, stabilizer, weights=[1, 1]):
    """
    Finds a matching to fix the errors in a 3D planar code given the positions.

    Parameters:
    -----------
    size -- The dimension of the code
    anyon -- A list of the locations of all '-1' value stabilizers in the 3D
             parity lattice. [[x0,y0,t0],[x1,y1,t1],...]
    stabilizer -- The stabilizer basis, can take the value "star" or "plaq"
    weights -- The multiplicative weighting that should be assigned to graph
               edges in the [space,time] dimensions. Default: [1,1]

    Returns:
    --------
    A list containing all the input anyon positions grouped into pairs.
    [[[x0,y0,t0],[x1,y1,t1]],[[x2,y2,t2],...

    """
    # NOTE: Use max time separation?
    # max_time_separation = 10

    N = len(anyons)
    if N == 0:
        return []

    # Append virtal anyons
    virtual_space = add_virtual_space(size, anyons, stabilizer)
    anyons = np.concatenate((anyons, virtual_space))

    nodes1, nodes2, weights = make_nodes_planar(anyons, weights)
    print(nodes1, nodes2, weights)
    if len(weights) == 0:
        return []
    matching = pm.getMatching_fast(2*N, nodes1, nodes2, weights)
    # REFORMAT MATCHING PAIRS
    # Take <matching> and turn it into a list of paired anyon positions.
    pairs_ind = [[i, matching[i]] for i in range(2*N) if matching[i] > i]

    pairs = [] if len(pairs_ind) == 0 else [anyons[p] for p in pairs_ind]

    pairs = np.array(pairs)
    # Remove unwanted pairs
    # Get the indices of the wanted pairs
    if stabilizer == "star":
        c = 1
    else:
        c = 0

    # Get the pairs where both are outside
    out = np.logical_or(pairs[:, :, c] == -1,
                        pairs[:, :, c] == 2*size - 1)
    # Invert to get the indices of the rest
    ind = np.invert(np.prod(out, 1).astype(bool))
    pairs = pairs[ind]

    return np.array(pairs)

def add_virtual_space(size, anyons, stabilizer):
    N = len(anyons)
    virtual = anyons.copy()
    if stabilizer == "star":
        virtual[:, 1] = -1
        virtual[:, 1][np.where(anyons[:, 1] > size)[0]] = (2*size - 1)
    else:
        virtual[:, 0] = -1
        virtual[:, 0][np.where(anyons[:, 0] > size)[0]] = (2*size - 1)
    return virtual
    # return np.concatenate((anyons, virtual))

def add_virtual_time(total_time, anyons):
    N = len(anyons)
    t = int(total_time/2)
    virtual = anyons.copy()
    virtual[:, 2] = -1
    virtual[:, 2][np.where(anyons[:, 2] > t)[0]] = time + 1

def make_nodes_planar(nodes, edges_weights=[1, 1]):
    nodes1 = []
    nodes2 = []
    weights = []

    N = int(len(nodes)/2)
    ws, wt = edges_weights
    wb = ws
    # Steps 1: Complete graph between all real nodes and real and virtual nodes
    for i in range(N):
        px, py, pt = nodes[i]

        for j in range(i+1, 2*N):
            qx, qy, qt = nodes[j]

            difft = (qt - pt)*wt
            diffx = abs(qx - px)*ws
            diffy = abs(qy - py)*ws
            weight = diffx + diffy + difft

            nodes1 += [i]
            nodes2 += [j]
            weights += [weight]

    for i in range(N, 2*N - 1):
        px, py, pt = nodes[i]

        for j in range(i+1, 2*N):
            qx, qy, qt = nodes[j]
            wt = (qt-pt)

            # if wt >= 10:
            #     break

            nodes1 += [i]
            nodes2 += [j]
            weights += [0]

    return nodes1, nodes2, weights
