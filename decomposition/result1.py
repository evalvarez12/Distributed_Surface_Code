"""
File to generate a analysis of how the error types depend on the
Fidelity of the GHZ state used into making the stabilizer measurements
in the distributed surface code.
"""
import numpy as np
import matplotlib.pyplot as plt
import protocols_det
import noise_modeling
import error_models as errs

# Initial parameters
ps = 0.0
pm = 0.015
pg = 0.015
pn = 0.1
system_size = 4
parity = "X"

# Initialize objects
model = noise_modeling.NoiseModel(system_size, parity)
model.separate_basis_parity()
prot = protocols_det.ProtocolsDeterministic(ps, pm, pg, pm)

# Choi state for noise noise modeling
choi = model._choi_state_ket(system_size)
choi = choi * choi.dag()
targets = list(range(system_size))

I_OK = []
I_NOK = []
E_OK = []
E_NOK = []

for f in np.arange(0, 1, .05):
    ghz = errs.generate_noisy_ghz(f, system_size)
    ps, rhos = prot.measure_ghz_stabilizer(choi, ghz, targets, parity)
    model.set_rho(rhos, ps)
    model.make_chi_matrix()

    I_OK += [model.chi["IIII_OK"]]
    I_NOK += [model.chi["IIII_NOK"]]
    E_OK += [model.chi["IXII_OK"]]
    E_NOK += [model.chi["IXII_NOK"]]

    # model.reset_chi()

plt.plot(I_OK)
plt.plot(I_NOK)
plt.plot(E_OK)
plt.plot(E_NOK)
plt.show()
