import noise_modeling
import protocols


# Set error values
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.1

# Initial parameters
system_size = 4
parity = "X"

# Initialize objects
prot = protocols.Protocols(ps, pm, pg, pn)
# Model takes the superoperator as a parameter
model = noise_modeling.NoiseModel(system_size, prot.monolithic, parity)

# Calculate chi matrix
model.separate_basis_parity()
model.make_chi_matrix()
