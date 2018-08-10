import ctypes
from numpy.ctypeslib import ndpointer
import os


def blossom_match(number_nodes, graph):

    number_edges = len(graph)

    PMlib = ctypes.CDLL("%s/PMlib.so" % "/".join((os.path.realpath(__file__)).split("/")[:-1]))
    PMlib.blossom_match.argtypes = [ctypes.c_int, ctypes.c_int,
                                 ctypes.POINTER(ctypes.c_int),
                                 ctypes.POINTER(ctypes.c_int),
                                 ctypes.POINTER(ctypes.c_int)]
    PMlib.blossom_match.restype = ndpointer(dtype=ctypes.c_int,
                                         shape=(number_nodes,))

    # initialize ctypes array and fill with edge data
    # nodes1 = graph[:, 0].ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    # nodes2 = graph[:, 1].ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    # weights = graph[:, 2].ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    nodes1 = (ctypes.c_int*number_edges)(*graph[:, 0])
    nodes2 = (ctypes.c_int*number_edges)(*graph[:, 1])
    weights = (ctypes.c_int*number_edges)(*graph[:, 2])

    # for i in range(number_edges):
    #     nodes1[i] = graph[i][0]
    #     nodes2[i] = graph[i][1]
    #     weights[i] = graph[i][2]
    #
    result = PMlib.blossom_match(ctypes.c_int(number_nodes),
                                 ctypes.c_int(number_edges),
                                 nodes1, nodes2, weights)

    return result
