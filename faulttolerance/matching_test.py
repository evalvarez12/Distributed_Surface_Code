"""
Test routines for matching functions

created on: 25/07/17

@author: eduardo
"""
import unittest
import numpy as np
import numpy.testing as nptest
import matching


class TestMatching(unittest.TestCase):

    def test_make_graph(self):
        nodes = np.random.randint(20, size=(100, 3))
        size = 10
        for i in range(100):
            graph = matching.make_graph_toric(size, nodes)


if __name__ == '__main__':
    unittest.main()
