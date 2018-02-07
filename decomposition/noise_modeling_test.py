import noise_modeling as nm
import protocols_tools as pt

# Initial parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.0
system_size = 4
parity = "X"

# Initialize objects
protocol = pt.ProtocolsDeterministic(ps, pm, pg, pn)
model = nm.NoiseModel(system_size, parity)

# Define function and apply superoperator
superoperator_function = protocol.local_stabilizer
model.apply_superoperator(superoperator_function)

# Separate basis using the parity symmetry
model.separate_basis_parity()

# Calculate chi matrix
model.make_chi_matrix()



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
