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

    if len(anyons) == 0:
        return []

    nodes1, nodes2, weights = make_nodes_planar(size, anyons,
                                                stabilizer,
                                                weights)
    N = len(weights)
    print("NUMBER nodes:", N)
    print("NODES:")
    print(nodes1, nodes2, weights)

    matching = pm.getMatching_fast(N, nodes1, nodes2, weights)
    # REFORMAT MATCHING PAIRS
    # Take <matching> and turn it into a list of paired anyon positions.
    pairs_ind = [[i, matching[i]] for i in range(N) if matching[i] > i]

    pairs = [] if len(pairs_ind) == 0 else [anyons[p] for p in pairs_ind]

    return pairs

def make_nodes_planar(size, nodes, stabilizer, edges_weights=[1, 1]):
    nodes1 = []
    nodes2 = []
    weights = []

    N = len(nodes)
    ws, wt = edges_weights
    wb = ws
    # Steps 1: Complete graph between all real nodes
    for i in range(N - 1):
        px, py, pt = nodes[i]

        for j in range(i+1, N):
            qx, qy, qt = nodes[j]

            difft = (qt - pt)*wt
            diffx = abs(qx - px)*ws
            diffy = abs(qy - py)*ws
            weight = diffx + diffy + difft

            nodes1 += [i]
            nodes2 += [j]
            weights += [weight]

    # Step 2: Generate list of boundary nodes linked to each real node

    boundary_nodes = []

    if stabilizer == "star":
        for i in range(N):
            px, py, pt = nodes[i]
            bx, bt = px, pt
            by = -1 if py < size else (2*size + 1)

            # difft = (qt - pt)*wt
            diffx = abs(bx - px)*ws
            diffy = abs(by - py)*ws
            weight = diffx + diffy

            nodes1 += [i]
            nodes2 += [i + N]
            weights += [weight]

            boundary_nodes += [(bx, by, bt)]

    elif stabilizer == "plaq":
        for i in range(N):
            px, py, pt = nodes[i]
            bx, bt = px, pt
            bx = -1 if px < size else (2*size + 1)

            # difft = (qt - pt)*wt
            diffx = abs(bx - px)*ws
            diffy = abs(by - py)*ws
            weight = diffx + diffy

            nodes1 += [i]
            nodes2 += [i + N]
            weights += [int(weight*wb/ws)]

            boundary_nodes += [(bx, by, bt)]

    # Step 3: Complete graph between all boundary nodes

    for i in range(N - 1):
        px, py, pt = boundary_nodes[i]

        for j in range(i+1, N):
            qx, qy, qt = boundary_nodes[j]
            wt = (qt-pt)

            if wt >= 5:
                break

            nodes1 += [N+i]
            nodes2 += [N+j]
            weights += [0]

    return nodes1, nodes2, weights
