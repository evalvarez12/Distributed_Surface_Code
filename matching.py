import blossom5.pyMatch as pm
import numpy as np


def match(size, anyons, surface, stabilizer, time, weights=[1, 1]):
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
    A list containing all the input anyon positions grouped into pairs.
    [[[x0,y0,t0],[x1,y1,t1]],[[x2,y2,t2],...
    """
    N_real = len(anyons)
    if N_real == 0:
        return []

    # Append virtal anyons
    cyclic = True
    if surface == "planar":
        anyons = add_virtual_space(size, anyons, stabilizer)
        cyclic = False
    if time != 0:
        anyons = add_virtual_time(time, anyons)
    N = len(anyons)

    graph = make_graph(size, anyons, N_real, weights, cyclic=cyclic)

    if N % 2 == 1:
        raise ValueError("Number of nodes is odd!")
    matching = pm.getMatching(N, graph)

    pairs_ind = np.array([[i, matching[i]] for i in range(N) if matching[i] > i])
    pairs = anyons[pairs_ind]
    # Pairs format:
    # np.array([[pair1, pair2], [pair1, pair2], ...])

    # Remove unwanted pairs
    if surface == "toric":
        pairs = pairs_remove_out_toric(time, pairs)
    elif surface == "planar":
        pairs = pairs_remove_out_planar(size, time, stabilizer, pairs)
    return pairs


def make_graph(size, nodes, N_real, weights=[1, 1], cyclic=True):
    # Nodes array format is:
    # [[x1, y1, t1], [x2, y2, t2], [x3, y3, t3], ...]
    N = len(nodes)
    # Spatial and time weights
    ws, wt = weights

    # NOTE: Cool method for making graph!!
    # Copy the nodes to obtain a matrix
    matrix = np.stack((nodes,)*N)
    # Get the substractions to find all distances
    matrix = matrix - np.expand_dims(nodes, 1)
    matrix = np.abs(matrix)
    # Use a mask to the correct distances
    # Get the upper triangular part to eliminate duplicates
    mask = np.triu(np.ones((N, N)))
    np.fill_diagonal(mask, 0)
    mask = mask.astype(bool)
    # Get the weights from the matrix of distances
    weights = matrix[mask]
    weights = np.abs(weights)
    # Calculate the min distance around the toroid and multiply by the weights
    if cyclic:
        # m is used to find shortest distance accros the toric
        m = 2*size
        weights[:, 0] = np.minimum(weights[:, 0], m - weights[:, 0])*ws
        weights[:, 1] = np.minimum(weights[:, 1], m - weights[:, 1])*ws
    else:
        weights[:, 0] *= ws
        weights[:, 1] *= ws

    weights[:, 2] *= wt
    # Add extra dimension to the weights to make concatenation possible
    weights = np.expand_dims(np.sum(weights, 1), 1)
    # The indices of the graph correspond to the indices of the matrix
    graph = np.array(np.where(mask)).transpose()
    # Join all the indexes and the weights
    graph = np.concatenate((graph, weights), 1)
    # Make weight between all virtual nodes = 0
    graph[graph[:, 0] >= N_real, 2] = 0

    return graph


def pairs_remove_out_planar(size, total_time, stabilizer, pairs):
    if stabilizer == "star":
        c = 1
    else:
        c = 0

    # Get the pairs where both are outside in space
    out = np.logical_or(pairs[:, :, c] == -1,
                        pairs[:, :, c] == 2*size - 1)
    # Invert to get the indices of the rest
    ind = np.invert(np.prod(out, 1).astype(bool))
    pairs = pairs[ind]

    # Get the pairs where both are outside in time
    out = np.prod(pairs[:, :, 2] == total_time + 1, 1)
    # Invert to get the indices of the rest
    ind = np.invert(out.astype(bool))
    pairs = pairs[ind]

    # Get the pairs where both are outside in space-time
    out = np.logical_or(pairs[:, :, 2] == total_time + 1,
                        pairs[:, :, c] == 2*size - 1)
    # Invert to get the indices of the rest
    ind = np.invert(np.prod(out, 1).astype(bool))
    pairs = pairs[ind]

    # Get the pairs where both are outside in space-time
    out = np.logical_or(pairs[:, :, c] == -1,
                        pairs[:, :, 2] == total_time + 1)
    # Invert to get the indices of the rest
    ind = np.invert(np.prod(out, 1).astype(bool))
    pairs = pairs[ind]
    return pairs


def pairs_remove_out_toric(total_time, pairs):
    # Get the pairs where both are outside in time
    out = np.prod(pairs[:, :, 2] == total_time + 1, 1)
    # Invert to get the indices of the rest
    ind = np.invert(out.astype(bool))
    pairs = pairs[ind]

    return pairs


def add_virtual_space(size, anyons, stabilizer):
    N = len(anyons)
    virtual = anyons.copy()
    if stabilizer == "star":
        virtual[:, 1] = -1
        virtual[:, 1][np.where(anyons[:, 1] > size)[0]] = (2*size - 1)
    else:
        virtual[:, 0] = -1
        virtual[:, 0][np.where(anyons[:, 0] > size)[0]] = (2*size - 1)
    return np.concatenate((anyons, virtual))


def add_virtual_time(total_time, anyons):
    N = len(anyons)
    t = int(total_time/2)
    virtual = anyons.copy()
    virtual[:, 2] = total_time + 1
    # virtual[:, 2][np.where(anyons[:, 2] > t)[0]] = total_time
    return np.concatenate((anyons, virtual))
