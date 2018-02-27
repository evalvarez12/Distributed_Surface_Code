import noise_modeling as nm
import stabilizer as st

# Initial parameters
ps = 0.006
pm = 0.006
pg = 0.006
system_size = 4
parity = "X"

# Initialize objects
stab = st.Stabilizer(ps, pm, pg)
model = nm.NoiseModel(system_size, parity)

# Define function and apply superoperator
superoperator_function = stab.local_stabilizer
model.apply_superoperator(superoperator_function)

# Separate basis using the parity symmetry
model.separate_basis_parity()

# Calculate chi matrix
model.make_chi_matrix()

print("Total sum: ", model.check_total_sum())
# Lenght must be 38
print("len: ", len(model.chi))


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
