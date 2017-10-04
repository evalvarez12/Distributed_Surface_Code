import unittest
import qutip as qt
import numpy as np
import numpy.testing as nptest
import pauli_basis as pb

# # Define variables for tests
# N = 3
# pos = 1


def assert_Equivalent(A, B):
    """
    Assert equivalence up to a global_phase between the two operators.
    """
    # Numpy array representations
    a = A.full()
    b = B.full()

    # Get the first non zero element
    non_zero_ind = np.nonzero(a)
    x = non_zero_ind[0][0]
    y = non_zero_ind[1][0]

    # Find the phase if theres any
    phase = b[x, y] / a[x, y]
    if phase * A == B:
        return True
    return False



class TestPauliBasis(unittest.TestCase):

    def test_string_to_operator(self):
        print("----------Test string_to_operator-----------")
        symbol = "IZYZ"
        op_ref = qt.sigmaz()
        op_ref = qt.tensor(qt.sigmay(), op_ref)
        op_ref = qt.tensor(qt.sigmaz(), op_ref)
        op_ref = qt.tensor(qt.qeye(2), op_ref)

        op = pb.string_to_operator(symbol, [0, 1, 2, 3], len(symbol))
        self.assertTrue(assert_Equivalent(op, op_ref))

    def test_symbol_product(self):
        print("----------Test symbol_product-----------")
        symA = "IXYZ"
        symB = "XYXZ"

        res = pb.symbol_product(symA, symB)
        self.assertEqual(res, "XZZI")

    def test_combinations(self):
        system_size = 4
        basis = pb.get_basis(list(range(system_size)), system_size)
        print(len(basis))


if __name__ == '__main__':
    unittest.main()
