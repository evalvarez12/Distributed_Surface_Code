import unittest
import qutip as qt
import numpy as np
import numpy.testing as nptest
import noise_modeling as nm

system_size = 2
superoperator_function = lambda x: x
model = nm.NoiseModel(system_size, superoperator_function)

choi = model.choi_state_ket(2)
print(choi)

symmetry = "Z"
print(model.pauli_basis.keys())
model.remove_sym_pauli_basis("X")
print(model.pauli_basis.keys())

# class TestNoiseModeling(unittest.TestCase):
#
#     def test_choi_state_ket(self):
#         print("----------Test choi_state_ket-----------")
#         choi = model.choi_state_ket(2)
#         print(choi)
#         # self.assertTrue(assert_Equivalent(op, op_ref))
#
#     def test_symbol_product(self):
#         print("----------Test symbol_product-----------")
#         symA = "IXYZ"
#         symB = "XYXZ"
#
#         res = pb.symbol_product(symA, symB)
#         self.assertEqual(res, "XZZI")
#
#
# if __name__ == '__main__':
#     unittest.main()
